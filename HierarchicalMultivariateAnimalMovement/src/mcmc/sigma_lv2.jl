

function sample_sigma_lv2_type2!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, sigmaPar::PosDefMatInverseWishart,parHier::OptionHierarchicalParameters,Clusterization::Vector{Clusterization_HDPHMM})

    kmax        = Likelihood[1].kmax

    Likelihood[1].clusterization.nonemptyC


    ###
    nanim_v2    = size(Likelihood)[1]
    nanim       = nanim_v2
    kmax        = Likelihood[1].kmax
    nt          = Likelihood[1].data.nt

    Aindex      = Matrix{Int16}(undef,kmax,5*nanim)
    
    for ian in 1:nanim_v2


        @inbounds Aindex[:,(ian-1)*5+1] = deepcopy(parHier.h_mu.clust[:,ian])
        @inbounds Aindex[:,(ian-1)*5+2] = deepcopy(parHier.h_eta.clust[:,ian])
        @inbounds Aindex[:,(ian-1)*5+3] = deepcopy(parHier.h_nu.clust[:,ian])
        @inbounds Aindex[:,(ian-1)*5+4] = deepcopy(parHier.h_rho.clust[:,ian])
        @inbounds Aindex[:,(ian-1)*5+5] = deepcopy(parHier.h_sigma.clust[:,ian])

    end

    for ian in 1:nanim_v2


        for k1 in 1:Likelihood[ian].clusterization.n_nonemptyC[1]
            @inbounds k = Likelihood[ian].clusterization.nonemptyC[k1]

            @inbounds k_beta  = parHier.h_sigma.clust[k,ian]
            @inbounds pp      = parHier.h_sigma.prob[k_beta]
            u       = rand(rng,Uniform(0.0,pp))

            # temporal indices

            @inbounds t_index = findall(x->x == k, Likelihood[ian].clusterization.zeta[1:(nt-1)])

            probvec = zeros(Float64,kmax)
            indexj  = zeros(Int16,kmax)
            n_j     = 0
            # uuu = rand(Uniform(0.0,1.0))
            # if uuu>0.95
            #     u = -0.1
            #     for j1 in 1:kmax
            #         @inbounds probvec[j1] += Float64(sum(Clusterization[ian].MatM[:,k])-1.0)*log(parHier.h_sigma.prob[j1])
            #     end
            # end
            for j1 in 1:kmax
                @inbounds Aindex[k,(ian-1)*5+5] = j1
                if (parHier.h_sigma.prob[j1]>u) & (check_unique(Aindex,k,kmax))

                    n_j         += one(n_j)
                    @inbounds indexj[n_j] = j1

                    @inbounds Likelihood[ian].sigma.parameteracc[k]    = deepcopy(parHier.h_sigma.par.parameteracc[j1])
                    @inbounds Likelihood[ian].sigma.parameteraccinv[k] = deepcopy(parHier.h_sigma.par.parameteraccinv[j1])

                    for t in t_index
                        @inbounds probvec[n_j] += compute_logpdf(Likelihood[ian], k, Int16(t+1))
                    end
                    #@inbounds probvec[n_j] += Float64(sum(Clusterization[ian].MatM[:,k])-1.0)*log(parHier.h_sigma.prob[j1])
                end
            end

            ## sample
            @inbounds sampvec = probvec[1:n_j]
            @inbounds samp_c = indexj[sample_discretevar(rng,sampvec)]

            @inbounds Likelihood[ian].sigma.parameteracc[k]    = deepcopy(parHier.h_sigma.par.parameteracc[samp_c])
            @inbounds Likelihood[ian].sigma.parameteraccinv[k] = deepcopy(parHier.h_sigma.par.parameteraccinv[samp_c])

            @inbounds Likelihood[ian].sigma.parameterprop[k]    = deepcopy(parHier.h_sigma.par.parameteracc[samp_c])
            @inbounds Likelihood[ian].sigma.parameterpropinv[k] = deepcopy(parHier.h_sigma.par.parameteraccinv[samp_c])
            @inbounds parHier.h_sigma.clust[k,ian]          = samp_c

            @inbounds Aindex[k,(ian-1)*5+5]                        = samp_c
        end

        if Likelihood[ian].clusterization.n_emptyC[1]>0
            for k1 in 1:Likelihood[ian].clusterization.n_emptyC[1]
                @inbounds k = Likelihood[ian].clusterization.emptyC[k1]

                probvec = zeros(Float64,kmax)
                indexj  = zeros(Int16,kmax)
                n_j     = 0
                for j1 in 1:kmax
                    @inbounds Aindex[k,(ian-1)*5+5] = j1
                    if check_unique(Aindex,k,kmax)
                        n_j          += one(n_j)
                        @inbounds indexj[n_j]  = j1
                        @inbounds probvec[n_j] = parHier.h_sigma.prob[j1]
                    end
                end
                @inbounds sampvec = probvec[1:n_j]
                @inbounds samp_c = indexj[sample_discretevar(rng,sampvec)]

                @inbounds Likelihood[ian].sigma.parameteracc[k] = deepcopy(parHier.h_sigma.par.parameteracc[samp_c])
                @inbounds Likelihood[ian].sigma.parameteraccinv[k] = deepcopy(parHier.h_sigma.par.parameteraccinv[samp_c])
                @inbounds parHier.h_sigma.clust[k,ian]          = samp_c

                @inbounds Likelihood[ian].sigma.parameterprop[k] = deepcopy(parHier.h_sigma.par.parameteracc[samp_c])
                @inbounds Likelihood[ian].sigma.parameterpropinv[k] = deepcopy(parHier.h_sigma.par.parameteraccinv[samp_c])

                @inbounds Aindex[k,(ian-1)*5+5]                        = samp_c
            end
        end

    end


    return nothing

end



function sample_sigma_lv2!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, sigmaPar::PosDefMatInverseWishart,parHier::OptionHierarchicalParameters,Clusterization::Vector{Clusterization_HDPHMM})

    kmax        = Likelihood[1].kmax

    Likelihood[1].clusterization.nonemptyC


    ###
    nanim_v2    = size(Likelihood)[1]
    kmax        = Likelihood[1].kmax
    nt          = Likelihood[1].data.nt



    for ian in 1:nanim_v2

        Aindex      = Matrix{Int16}(undef,kmax,5)
        @inbounds Aindex[:,1] = deepcopy(parHier.h_mu.clust[:,ian])
        @inbounds Aindex[:,2] = deepcopy(parHier.h_eta.clust[:,ian])
        @inbounds Aindex[:,3] = deepcopy(parHier.h_nu.clust[:,ian])
        @inbounds Aindex[:,4] = deepcopy(parHier.h_rho.clust[:,ian])
        @inbounds Aindex[:,5] = deepcopy(parHier.h_sigma.clust[:,ian])

        for k1 in 1:Likelihood[ian].clusterization.n_nonemptyC[1]
            @inbounds k = Likelihood[ian].clusterization.nonemptyC[k1]

            @inbounds k_beta  = parHier.h_sigma.clust[k,ian]
            @inbounds pp      = parHier.h_sigma.prob[k_beta]
            u       = rand(rng,Uniform(0.0,pp))

            # temporal indices

            @inbounds t_index = findall(x->x == k, Likelihood[ian].clusterization.zeta[1:(nt-1)])

            probvec = zeros(Float64,kmax)
            indexj  = zeros(Int16,kmax)
            n_j     = 0
            # uuu = rand(Uniform(0.0,1.0))
            # if uuu>0.95
            #     u = -0.1
            #     for j1 in 1:kmax
            #         @inbounds probvec[j1] += Float64(sum(Clusterization[ian].MatM[:,k])-1.0)*log(parHier.h_sigma.prob[j1])
            #     end
            # end
            for j1 in 1:kmax
                @inbounds Aindex[k,5] = j1
                if (parHier.h_sigma.prob[j1]>u) & (check_unique(Aindex,k,kmax))

                    n_j         += one(n_j)
                    @inbounds indexj[n_j] = j1

                    @inbounds Likelihood[ian].sigma.parameteracc[k]    = deepcopy(parHier.h_sigma.par.parameteracc[j1])
                    @inbounds Likelihood[ian].sigma.parameteraccinv[k] = deepcopy(parHier.h_sigma.par.parameteraccinv[j1])

                    for t in t_index
                        @inbounds probvec[n_j] += compute_logpdf(Likelihood[ian], k, Int16(t+1))
                    end
                    #@inbounds probvec[n_j] += Float64(sum(Clusterization[ian].MatM[:,k])-1.0)*log(parHier.h_sigma.prob[j1])
                end
            end

            ## sample
            @inbounds sampvec = probvec[1:n_j]
            @inbounds samp_c = indexj[sample_discretevar(rng,sampvec)]

            @inbounds Likelihood[ian].sigma.parameteracc[k]    = deepcopy(parHier.h_sigma.par.parameteracc[samp_c])
            @inbounds Likelihood[ian].sigma.parameteraccinv[k] = deepcopy(parHier.h_sigma.par.parameteraccinv[samp_c])

            @inbounds Likelihood[ian].sigma.parameterprop[k]    = deepcopy(parHier.h_sigma.par.parameteracc[samp_c])
            @inbounds Likelihood[ian].sigma.parameterpropinv[k] = deepcopy(parHier.h_sigma.par.parameteraccinv[samp_c])
            @inbounds parHier.h_sigma.clust[k,ian]          = samp_c

            @inbounds Aindex[k,5]                        = samp_c
        end

        if Likelihood[ian].clusterization.n_emptyC[1]>0
            for k1 in 1:Likelihood[ian].clusterization.n_emptyC[1]
                @inbounds k = Likelihood[ian].clusterization.emptyC[k1]

                probvec = zeros(Float64,kmax)
                indexj  = zeros(Int16,kmax)
                n_j     = 0
                for j1 in 1:kmax
                    @inbounds Aindex[k,5] = j1
                    if check_unique(Aindex,k,kmax)
                        n_j          += one(n_j)
                        @inbounds indexj[n_j]  = j1
                        @inbounds probvec[n_j] = parHier.h_sigma.prob[j1]
                    end
                end
                @inbounds sampvec = probvec[1:n_j]
                @inbounds samp_c = indexj[sample_discretevar(rng,sampvec)]

                @inbounds Likelihood[ian].sigma.parameteracc[k] = deepcopy(parHier.h_sigma.par.parameteracc[samp_c])
                @inbounds Likelihood[ian].sigma.parameteraccinv[k] = deepcopy(parHier.h_sigma.par.parameteraccinv[samp_c])
                @inbounds parHier.h_sigma.clust[k,ian]          = samp_c

                @inbounds Likelihood[ian].sigma.parameterprop[k] = deepcopy(parHier.h_sigma.par.parameteracc[samp_c])
                @inbounds Likelihood[ian].sigma.parameterpropinv[k] = deepcopy(parHier.h_sigma.par.parameteraccinv[samp_c])

                @inbounds Aindex[k,5]                        = samp_c
            end
        end

    end


    return nothing

end



function sample_sigma_split_merge!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, sigmaPar::PosDefMatInverseWishart,parHier::OptionHierarchicalParameters,Clusterization::Vector{Clusterization_HDPHMM})

    kmax        = Likelihood[1].kmax


    #print("start" , "\n")
    ###
    nanim_v2    = size(Likelihood)[1]
    kmax        = Likelihood[1].kmax
    nt          = Likelihood[1].data.nt

    # for k in 1:kmax
    #     parHier.h_sigma.par.parameteraccinv[k] = inv(parHier.h_sigma.par.parameteracc[k])
    # end
    nsamp       = Int16(nanim_v2)

    ind_a1      = zeros(Int16,nsamp)
    ind_a2      = zeros(Int16,nsamp)
    ind_k1      = zeros(Int16,nsamp)
    ind_k2      = zeros(Int16,nsamp)
    ind_j1      = zeros(Int16,nsamp)
    ind_j2      = zeros(Int16,nsamp)

    for is in 1:nsamp

        ind_a1[is] = Int16(rand(1:nanim_v2,1)[1])
        ind_a2[is] = Int16(rand(1:nanim_v2,1)[1])

        ind_k1[is] = Int16(rand(1:Likelihood[ind_a1[is]].clusterization.n_nonemptyC[1],1)[1])
        ind_k2[is] = Int16(rand(1:Likelihood[ind_a2[is]].clusterization.n_nonemptyC[1],1)[1])

        ind_k1[is] = Likelihood[ind_a1[is]].clusterization.nonemptyC[ind_k1[is]]
        ind_k2[is] = Likelihood[ind_a2[is]].clusterization.nonemptyC[ind_k2[is]]

        ind_j1[is] = parHier.h_sigma.clust[ind_k1[is],ind_a1[is]]
        ind_j2[is] = parHier.h_sigma.clust[ind_k2[is],ind_a2[is]]

    end

    for is in 1:nsamp


        if ind_j1[is]!=ind_j2[is]

            probvec   = zeros(Float64,4)

            t_index_1 = findall(x->x == ind_k1[is], Likelihood[ind_a1[is]].clusterization.zeta[1:(nt-1)])
            t_index_2 = findall(x->x == ind_k2[is], Likelihood[ind_a2[is]].clusterization.zeta[1:(nt-1)])

            prob1       = Float64(0.0)
            prob2       = Float64(0.0)
            prob1_new   = Float64(0.0)
            prob2_new   = Float64(0.0)


            Likelihood[ind_a1[is]].sigma.parameteracc[ind_k1[is]] = deepcopy(parHier.h_sigma.par.parameteracc[ind_j1[is]])
            Likelihood[ind_a1[is]].sigma.parameteraccinv[ind_k1[is]] = deepcopy(parHier.h_sigma.par.parameteraccinv[ind_j1[is]])
            n1 = Int16(0)
            for t in t_index_1
                @inbounds prob1 += compute_logpdf(Likelihood[ind_a1[is]], ind_k1[is], Int16(t+1))
                n1 += one(n1)
            end

            Likelihood[ind_a2[is]].sigma.parameteracc[ind_k2[is]] = deepcopy(parHier.h_sigma.par.parameteracc[ind_j2[is]])
            Likelihood[ind_a2[is]].sigma.parameteraccinv[ind_k2[is]] = deepcopy(parHier.h_sigma.par.parameteraccinv[ind_j2[is]])
            n2 = Int16(0)
            for t in t_index_2
                @inbounds prob2 += compute_logpdf(Likelihood[ind_a2[is]], ind_k2[is], Int16(t+1))
                n2 += one(n2)
            end

            Likelihood[ind_a1[is]].sigma.parameteracc[ind_k1[is]] = deepcopy(parHier.h_sigma.par.parameteracc[ind_j2[is]])
            Likelihood[ind_a1[is]].sigma.parameteraccinv[ind_k1[is]] = deepcopy(parHier.h_sigma.par.parameteraccinv[ind_j2[is]])
            for t in t_index_1
                @inbounds prob1_new += compute_logpdf(Likelihood[ind_a1[is]], ind_k1[is], Int16(t+1))
            end

            Likelihood[ind_a2[is]].sigma.parameteracc[ind_k2[is]] = deepcopy(parHier.h_sigma.par.parameteracc[ind_j1[is]])
            Likelihood[ind_a2[is]].sigma.parameteraccinv[ind_k2[is]] = deepcopy(parHier.h_sigma.par.parameteraccinv[ind_j1[is]])
            for t in t_index_2
                @inbounds prob2_new += compute_logpdf(Likelihood[ind_a2[is]], ind_k2[is], Int16(t+1))
            end

            probvec[1] = prob1+prob2
            probvec[1] += Float64(n1)*log(parHier.h_sigma.prob[ind_j1[is]])
            probvec[1] += Float64(n2)*log(parHier.h_sigma.prob[ind_j2[is]])

            probvec[2] = prob1_new+prob2_new
            probvec[2] += Float64(n1)*log(parHier.h_sigma.prob[ind_j2[is]])
            probvec[2] += Float64(n2)*log(parHier.h_sigma.prob[ind_j1[is]])

            probvec[3] = prob1+prob2_new
            probvec[3] += Float64(n1)*log(parHier.h_sigma.prob[ind_j1[is]])
            probvec[3] += Float64(n2)*log(parHier.h_sigma.prob[ind_j1[is]])

            probvec[4] = prob1_new+prob2
            probvec[4] += Float64(n1)*log(parHier.h_sigma.prob[ind_j2[is]])
            probvec[4] += Float64(n2)*log(parHier.h_sigma.prob[ind_j2[is]])


            samp_c = sample_discretevar(rng,probvec)
            #print("sample ", samp_c, "\n")
            #print(probvec[1]," ",probvec[2]," ",probvec[3]," ",probvec[4],"\n")
            if samp_c==1

                Likelihood[ind_a1[is]].sigma.parameteracc[ind_k1[is]]    = deepcopy(parHier.h_sigma.par.parameteracc[ind_j1[is]])
                Likelihood[ind_a1[is]].sigma.parameteraccinv[ind_k1[is]] = deepcopy(parHier.h_sigma.par.parameteraccinv[ind_j1[is]])
                parHier.h_sigma.clust[ind_k1[is],ind_a1[is]]          = ind_j1[is]

                Likelihood[ind_a2[is]].sigma.parameteracc[ind_k2[is]]    = deepcopy(parHier.h_sigma.par.parameteracc[ind_j2[is]])
                Likelihood[ind_a2[is]].sigma.parameteraccinv[ind_k2[is]] = deepcopy(parHier.h_sigma.par.parameteraccinv[ind_j2[is]])
                parHier.h_sigma.clust[ind_k2[is],ind_a2[is]]          = ind_j2[is]

            elseif samp_c==2

                Likelihood[ind_a1[is]].sigma.parameteracc[ind_k1[is]]    = deepcopy(parHier.h_sigma.par.parameteracc[ind_j2[is]])
                Likelihood[ind_a1[is]].sigma.parameteraccinv[ind_k1[is]] = deepcopy(parHier.h_sigma.par.parameteraccinv[ind_j2[is]])
                parHier.h_sigma.clust[ind_k1[is],ind_a1[is]]          = ind_j2[is]

                Likelihood[ind_a2[is]].sigma.parameteracc[ind_k2[is]]    = deepcopy(parHier.h_sigma.par.parameteracc[ind_j1[is]])
                Likelihood[ind_a2[is]].sigma.parameteraccinv[ind_k2[is]] = deepcopy(parHier.h_sigma.par.parameteraccinv[ind_j1[is]])
                parHier.h_sigma.clust[ind_k2[is],ind_a2[is]]          = ind_j1[is]

            elseif samp_c==3

                Likelihood[ind_a1[is]].sigma.parameteracc[ind_k1[is]]    = deepcopy(parHier.h_sigma.par.parameteracc[ind_j1[is]])
                Likelihood[ind_a1[is]].sigma.parameteraccinv[ind_k1[is]] = deepcopy(parHier.h_sigma.par.parameteraccinv[ind_j1[is]])
                parHier.h_sigma.clust[ind_k1[is],ind_a1[is]]          = ind_j1[is]

                Likelihood[ind_a2[is]].sigma.parameteracc[ind_k2[is]]    = deepcopy(parHier.h_sigma.par.parameteracc[ind_j1[is]])
                Likelihood[ind_a2[is]].sigma.parameteraccinv[ind_k2[is]] = deepcopy(parHier.h_sigma.par.parameteraccinv[ind_j1[is]])
                parHier.h_sigma.clust[ind_k2[is],ind_a2[is]]          = ind_j1[is]

            else
                Likelihood[ind_a1[is]].sigma.parameteracc[ind_k1[is]]    = deepcopy(parHier.h_sigma.par.parameteracc[ind_j2[is]])
                Likelihood[ind_a1[is]].sigma.parameteraccinv[ind_k1[is]] = deepcopy(parHier.h_sigma.par.parameteraccinv[ind_j2[is]])
                parHier.h_sigma.clust[ind_k1[is],ind_a1[is]]          = ind_j2[is]

                Likelihood[ind_a2[is]].sigma.parameteracc[ind_k2[is]]    = deepcopy(parHier.h_sigma.par.parameteracc[ind_j2[is]])
                Likelihood[ind_a2[is]].sigma.parameteraccinv[ind_k2[is]] = deepcopy(parHier.h_sigma.par.parameteraccinv[ind_j2[is]])
                parHier.h_sigma.clust[ind_k2[is],ind_a2[is]]          = ind_j2[is]
            end

        end
    end

    return nothing

end
