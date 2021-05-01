#### #### #### #### #### #### #### #### #### ####
#### #### #### sample mu0
#### #### #### #### #### #### #### #### #### ####

function sample_rho!(rng::MersenneTwister,Likelihood::AbstractLikelihood, rhoPar::AbstractVecPar)
    error(string("sample_rho not defined for Likelihood ", typeof(Likelihood), " and rho ", typeof(rhoPar)) )
end
function sample_rho!(rng::MersenneTwister,Likelihood::Likelihood_OU_CircLinmodel, rhoPar::VecParMvNoPrior)

end

function sample_rho!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, rhoPar::VecParMvUniform)

    nanim = size(Likelihood)[1]
    for ian in 1:nanim
        sample_rho!(rng,Likelihood[ian], Likelihood[ian].rho)
    end

end


function sample_rho!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, rhoPar::VecParMvUniform, parHier::OptionHierarchicalParameters)

    #Likelihood  = MCMCLikelihood
    # rhoPar      = MCMCLikelihood.rho

    nanim_v2    = size(Likelihood)[1]
    kmax        = Likelihood[1].kmax

    Loglike = zeros(Float64,kmax)


    prop_lv1 = [zeros(kmax) for _ in 1:nanim_v2]
    prop_lv2 = zeros(Float64,kmax)
    #lim1 = rhoPar.prior.v[1].a
    #lim2 = rhoPar.prior.v[1].b

    lim1 = 0.0
    lim2 = 1.0
    #error("controllare campionamento rho come l'altro articolo")
    prob_lim = Float64(1.0/3.0)*3.0


    VecK = zeros(Int16,kmax)
    for ian in 1:nanim_v2
        for k1 in 1:Likelihood[ian].clusterization.n_nonemptyC[1]
            k = Likelihood[ian].clusterization.nonemptyC[k1]

            VecK[parHier.h_rho.clust[k,ian]] = 1
        end
    end

    for k in 1:kmax
        if VecK[k] == 0
            uv   = rand(rng,Uniform(0.0,3.0))
            if uv < 1.0
                prop_lv2[k] = 0.0
            elseif uv > 2.0
                prop_lv2[k] = 1.0
            else
                prop_lv2[k] = rand(rng,Uniform(0.0,1.0))
            end

            if (prop_lv2[k]>0.0) && (prop_lv2[k]<1.0)
                Loglike[k] -=  log(1.0/3.0)
            else
                Loglike[k] -=  log(2.0/3.0)
            end

            if (parHier.h_rho.par.parameteracc[k][1]>0.0) && (parHier.h_rho.par.parameteracc[k][1]<1.0)
                Loglike[k] +=  log(1.0/3.0)
            else
                Loglike[k] +=  log(2.0/3.0)
            end

        else
            isinternal  = Int16(0)
            uv          = rand(rng,Uniform(0.0,3.0))
            if uv < prob_lim
                prop_lv2[k] = 0.0
            elseif uv > 3.0-prob_lim
                prop_lv2[k] = 1.0
            else
                u    = rand(rng,Uniform(0.0,1.0))
                if u<0.25
                    prop_lv2[k] = rand(Uniform(parHier.h_rho.par.parameteracc[k][1]-0.1, parHier.h_rho.par.parameteracc[k][1]+0.1))
                    #prop_lv2[k] = rand(rng,Normal(parHier.h_rho.par.parameteracc[k][1]-lim1, 0.1))
                elseif u<0.5
                    prop_lv2[k] = rand(Uniform(parHier.h_rho.par.parameteracc[k][1]-0.01, parHier.h_rho.par.parameteracc[k][1]+0.01))
                elseif u<0.75
                    prop_lv2[k] = rand(Uniform(parHier.h_rho.par.parameteracc[k][1]-0.001, parHier.h_rho.par.parameteracc[k][1]+0.001))
                else
                    prop_lv2[k] = rand(Uniform(parHier.h_rho.par.parameteracc[k][1]-0.0001, parHier.h_rho.par.parameteracc[k][1]+0.0001))
                end
                isinternal = 1
                #prop_lv2[k] = lim1+rem(prop_lv2[k],lim2-lim1 ,RoundDown)
            end

            if (prop_lv2[k]>0.0) && (prop_lv2[k]<1.0)
                Loglike[k] -= log(1.0-2.0*prob_lim/3.0)
                appP = 0.0
                if abs(prop_lv2[k]-parHier.h_rho.par.parameteracc[k][1])<0.1
                    appP += 0.25*1.0/(0.1*2.0)
                end
                if abs(prop_lv2[k]-parHier.h_rho.par.parameteracc[k][1])<0.01
                    appP += 0.25*1.0/(0.01*2.0)
                end
                if abs(prop_lv2[k]-parHier.h_rho.par.parameteracc[k][1])<0.001
                    appP += 0.25*1.0/(0.001*2.0)
                end
                if abs(prop_lv2[k]-parHier.h_rho.par.parameteracc[k][1])<0.0001
                    appP += 0.25*1.0/(0.0001*2.0)
                end
                Loglike[k] -= log(appP)
            else
                Loglike[k] -= log(prob_lim/3.0)
            end

            if (parHier.h_rho.par.parameteracc[k][1]>0.0) && (parHier.h_rho.par.parameteracc[k][1]<1.0)

                Loglike[k] += log(1.0-2.0*prob_lim/3.0)

                appP = 0.0
                if abs(prop_lv2[k]-parHier.h_rho.par.parameteracc[k][1])<0.1
                    appP += 0.25*1.0/(0.1*2.0)
                end
                if abs(prop_lv2[k]-parHier.h_rho.par.parameteracc[k][1])<0.01
                    appP += 0.25*1.0/(0.01*2.0)
                end
                if abs(prop_lv2[k]-parHier.h_rho.par.parameteracc[k][1])<0.001
                    appP += 0.25*1.0/(0.001*2.0)
                end
                if abs(prop_lv2[k]-parHier.h_rho.par.parameteracc[k][1])<0.0001
                    appP += 0.25*1.0/(0.0001*2.0)
                end

                Loglike[k] += log(appP)
            else
                Loglike[k] += log(prob_lim/3.0)
            end

            if isinternal == 1
                prop_lv2[k] = lim1+rem(prop_lv2[k],lim2-lim1 ,RoundDown)
                #prop[k] = lim1+rem(prop[k],lim2-lim1 ,RoundDown)
            end

        end
    end


    for ian in 1:nanim_v2
        for k in 1:kmax
            prop_lv1[ian][k] = prop_lv2[parHier.h_rho.clust[k,ian]]
        end
    end


    for ian in 1:nanim_v2



        nc          = Likelihood[ian].data.ncol
        nanim       = Likelihood[ian].data.nanimals
        nt          = Likelihood[ian].data.nt

        psi         = Likelihood[ian].nu.parameteracc
        rho         = Likelihood[ian].rho.parameteracc
        muC         = Likelihood[ian].eta.parameteracc
        mu0         = Likelihood[ian].mu.parameteracc
        sigmainv    = Likelihood[ian].sigma.parameteraccinv
        sigma       = Likelihood[ian].sigma.parameteracc

        zeta        = Likelihood[ian].clusterization.zeta
        Obs         = Likelihood[ian].data.data

        StartAngle  = Likelihood[ian].Angle.parameteracc




        prop        = prop_lv1[ian]
        for ianim in 1:nanim

            Dpsirho =  Vector{Matrix{Float64}}()
            for k in 1:kmax
                app  = repeat([1:nanim;], inner=2, outer=1)
                push!(Dpsirho, diagm(  psi[k][ app ].*(1.0 .-rho[k][ app ])   )   )
            end
            Drho =  Vector{Matrix{Float64}}()
            for k in 1:kmax
                app  = repeat([1:nanim;], inner=2, outer=1)
                push!(Drho, deepcopy(diagm(rho[k][ app ]  )))
            end



            MatR = zeros(Float64,nc,nc)
            MatR_prop = zeros(Float64,nc,nc)
            MatR2 = zeros(Float64,nc,nc)
            MatR2_prop = zeros(Float64,nc,nc)




            # for k in 1:kmax
            #     prop[k] = rho[k][ianim]+0.00001
            # end
            Dpsirho_prop     = deepcopy(Dpsirho)
            Drho_prop        = deepcopy(Drho)

            for k in 1:kmax
                WW = (ianim-1)*2 .+ [1;2]
                Dpsirho_prop[k][WW,WW]     = (1.0-prop[k]).*psi[k][ianim].*Matrix(I,2,2)
                Drho_prop[k][WW,WW]        = prop[k].*Matrix(I,2,2)
            end


            Cangle      = deepcopy(StartAngle[1])
            MatR        = zeros(Float64,nc,nc)
            MatRprop    = zeros(Float64,nc,nc)
            MatR2        = zeros(Float64,nc,nc)
            MatR2prop    = zeros(Float64,nc,nc)
            for i in 2:nt
                k       = zeta[i-1]
                yt1P    = Obs[i]  #view(Likelihood.data.data,2,:)
                yt      = Obs[i-1]

                for ianim2 in 1:nanim
                    W = [1,2] .+(ianim2-1)*2
                    MatR[W,W]       = [ cos(Cangle[ianim2]*rho[k][ianim2]) -sin(Cangle[ianim2]*rho[k][ianim2]); sin(Cangle[ianim2]*rho[k][ianim2])  cos(Cangle[ianim2]*rho[k][ianim2])  ]
                    MatR2[W,W]       = [ cos(Cangle[ianim2] ) -sin(Cangle[ianim2] ); sin(Cangle[ianim2] )  cos(Cangle[ianim2] )  ]
                    if ianim2!=ianim
                        MatR_prop[W,W]  = [ cos(Cangle[ianim2]*rho[k][ianim2]) -sin(Cangle[ianim2]*rho[k][ianim2]); sin(Cangle[ianim2]*rho[k][ianim2])  cos(Cangle[ianim2]*rho[k][ianim2])  ]
                        MatR2_prop[W,W]  = [ cos(Cangle[ianim2] ) -sin(Cangle[ianim2] ); sin(Cangle[ianim2] )  cos(Cangle[ianim2] )  ]
                    else
                        MatR_prop[W,W]  = [ cos(Cangle[ianim2]*prop[k]) -sin(Cangle[ianim2]*prop[k]); sin(Cangle[ianim2]*prop[k])  cos(Cangle[ianim2]*prop[k])  ]
                        MatR2_prop[W,W]  = [ cos(Cangle[ianim2]) -sin(Cangle[ianim2] ); sin(Cangle[ianim2] )  cos(Cangle[ianim2] )  ]
                    end

                end

                ##
                Mean_prop   = transpose(MatR_prop)*(yt+Dpsirho_prop[k]*(mu0[k]-yt) + Drho_prop[k]*MatR2_prop*muC[k])
                Mean_acc    = transpose(MatR)*(yt+Dpsirho[k]*(mu0[k]-yt) + Drho[k]*MatR2*muC[k])


                Loglike[parHier.h_rho.clust[k,ian]] += logpdf(MvNormal(Mean_prop, sigma[k]), transpose(MatR_prop)*yt1P)
                Loglike[parHier.h_rho.clust[k,ian]] -= logpdf(MvNormal(Mean_acc,  sigma[k]),   transpose(MatR)*yt1P)
                ##

                for j in 1:nanim
                    WW = [1 2] .+ (j-1)*2
                    Cangle[j] = atan(yt1P[WW][2]-yt[WW][2],yt1P[WW][1]-yt[WW][1])
                end
                #print(Loglike)
                #print("\n\n")

            end



            # print("LogLike")
            # print(Loglike)
            # print("\n")
        end

    end


    ianim = 1
    for k in 1:kmax
        u = rand(rng,Uniform(0.0,1.0))
        if u<exp(Loglike[k])
            parHier.h_rho.par.parameteracc[k][ianim]  = prop_lv2[k]
            parHier.h_rho.par.parameterprop[k][ianim] = prop_lv2[k]
        end
    end

    for ian in 1:nanim_v2
        rho         = Likelihood[ian].rho
        for k in 1:kmax
            rho.parameteracc[k][1]  = parHier.h_rho.par.parameteracc[parHier.h_rho.clust[k,ian]][1]
            rho.parameterprop[k][1]    = rho.parameteracc[k][1]
        end
    end


    return nothing

end


function sample_rho!(rng::MersenneTwister,Likelihood::Likelihood_OU_CircLinmodel, rhoPar::VecParMvUniform)

    #Likelihood  = MCMCLikelihood
    # rhoPar      = MCMCLikelihood.rho

    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    psi         = Likelihood.nu.parameteracc
    rho         = Likelihood.rho.parameteracc
    muC         = Likelihood.eta.parameteracc
    mu0         = Likelihood.mu.parameteracc
    sigmainv    = Likelihood.sigma.parameteraccinv
    sigma       = Likelihood.sigma.parameteracc

    zeta        = Likelihood.clusterization.zeta
    Obs         = Likelihood.data.data

    StartAngle  = Likelihood.Angle.parameteracc

    # CovInv_AgivenB      = Vector{Matrix{Float64}}()
    # for k in 1:kmax
    #     push!(CovInv_AgivenB, zeros(Float64,2,2))
    # end
    #
    # CovABInvCovB      = Vector{Matrix{Float64}}()
    # for k in 1:kmax
    #     push!(CovABInvCovB, zeros(Float64,2,nc-2))
    # end


    for ianim in 1:nanim

        Dpsirho =  Vector{Matrix{Float64}}()
        for k in 1:kmax
            app  = repeat([1:nanim;], inner=2, outer=1)
            push!(Dpsirho, diagm(  psi[k][ app ].*(1.0 .-rho[k][ app ])   )   )
        end
        Drho =  Vector{Matrix{Float64}}()
        for k in 1:kmax
            app  = repeat([1:nanim;], inner=2, outer=1)
            push!(Drho, deepcopy(diagm(rho[k][ app ]  )))
        end



        MatR = zeros(Float64,nc,nc)
        MatR_prop = zeros(Float64,nc,nc)

        MatR2 = zeros(Float64,nc,nc)
        MatR2_prop = zeros(Float64,nc,nc)

        Loglike = zeros(Float64,kmax)
        prop = zeros(Float64,kmax)
        lim1 = rhoPar.prior.v[ianim].a
        lim2 = rhoPar.prior.v[ianim].b



        # for k in 1:kmax
        #     uv   = rand(rng,Uniform(0.0,3.0))
        #     if uv < 1.0
        #         prop[k] = 0.0
        #     elseif uv > 2.0
        #         prop[k] = 1.0
        #     else
        #         u    = rand(rng,Uniform(0.0,1.0))
        #         if u<0.25
        #             prop[k] = rand(rng,Normal(rho[k][ianim]-lim1, 1.0))
        #         elseif u<0.5
        #             prop[k] = rand(rng,Normal(rho[k][ianim]-lim1, 0.1))
        #         elseif u<0.75
        #             prop[k] = rand(rng,Normal(rho[k][ianim]-lim1, 0.001))
        #         else
        #             prop[k] = rand(rng,Normal(rho[k][ianim]-lim1, 0.0001))
        #         end
        #
        #         prop[k] = lim1+rem(prop[k],lim2-lim1 ,RoundDown)
        #     end
        #
        #     #prop[k] = rho[k][ianim]
        # end
        prob_lim = Float64(1.0/3.0)*3.0
        for k in 1:kmax
            isinternal = Int16(0)
            uv   = rand(Uniform(0.0,3.0))
            if uv < prob_lim
                prop[k] = 0.0
            elseif uv > 3.0-prob_lim
                prop[k] = 1.0
            else
                u    = rand(Uniform(0.0,1.0))
                if u<0.25
                    prop[k] = rand(Uniform(rho[k][ianim]-0.5, rho[k][ianim]+0.5))
                elseif u<0.5
                    prop[k] = rand(Uniform(rho[k][ianim]-0.25, rho[k][ianim]+0.25))
                elseif u<0.75
                    prop[k] = rand(Uniform(rho[k][ianim]-0.1, rho[k][ianim]+0.1))
                else
                    prop[k] = rand(Uniform(rho[k][ianim]-0.01, rho[k][ianim]+0.01))
                end
                isinternal = 1
                #prop[k] = lim1+rem(prop[k],lim2-lim1 ,RoundDown)
            end


            if (prop[k]>0.0) && (prop[k]<1.0)
                Loglike[k] -= log(1.0-2.0*prob_lim/3.0)
                appP = 0.0
                if abs(prop[k]-rho[k][ianim])<0.5
                    appP += 0.25*1.0/(0.5*2.0)
                end
                if abs(prop[k]-rho[k][ianim])<0.25
                    appP += 0.25*1.0/(0.25*2.0)
                end
                if abs(prop[k]-rho[k][ianim])<0.1
                    appP += 0.25*1.0/(0.1*2.0)
                end
                if abs(prop[k]-rho[k][ianim])<0.01
                    appP += 0.25*1.0/(0.01*2.0)
                end
                Loglike[k] -= log(appP)
            else
                Loglike[k] -= log(prob_lim/3.0)
            end

            if (rho[k][ianim]>0.0) && (rho[k][ianim]<1.0)

                Loglike[k] += log(1.0-2.0*prob_lim/3.0)

                appP = 0.0
                if abs(prop[k]-rho[k][ianim])<0.5
                    appP += 0.25*1.0/(0.5*2.0)
                end
                if abs(prop[k]-rho[k][ianim])<0.25
                    appP += 0.25*1.0/(0.25*2.0)
                end
                if abs(prop[k]-rho[k][ianim])<0.1
                    appP += 0.25*1.0/(0.1*2.0)
                end
                if abs(prop[k]-rho[k][ianim])<0.01
                    appP += 0.25*1.0/(0.01*2.0)
                end

                Loglike[k] += log(appP)
            else
                Loglike[k] += log(prob_lim/3.0)
            end

            if isinternal == 1
                prop[k] = lim1+rem(prop[k],lim2-lim1 ,RoundDown)
            end


            #prop[k] = rho[k][ianim]
        end

        # prob_lim = Float64(0.33)*3
        # for k in 1:kmax
        #     uv   = rand(Uniform(0.0,3.0))
        #     if uv < prob_lim
        #         prop[k] = 0.0
        #     elseif uv > 3.0-prob_lim
        #         prop[k] = 1.0
        #     else
        #         u    = rand(Uniform(0.0,1.0))
        #         if u<0.25
        #             prop[k] = rand(Normal(rho[k][ianim]-lim1, 0.1))
        #         elseif u<0.5
        #             prop[k] = rand(Normal(rho[k][ianim]-lim1, 0.01))
        #         elseif u<0.75
        #             prop[k] = rand(Normal(rho[k][ianim]-lim1, 0.001))
        #         else
        #             prop[k] = rand(Normal(rho[k][ianim]-lim1, 0.0001))
        #         end
        #
        #         prop[k] = lim1+rem(prop[k],lim2-lim1 ,RoundDown)
        #     end
        #
        #
        #     if (prop[k]>0.0) && (prop[k]<1.0)
        #         Loglike[k] -= log(1.0-2.0*prob_lim/3.0)+ log(
        #         0.25*exp(logpdf(Normal(prop[k],0.1),rho[k][ianim]))+
        #         0.25*exp(logpdf(Normal(prop[k],0.01),rho[k][ianim]))+
        #         0.25*exp(logpdf(Normal(prop[k],0.001),rho[k][ianim]))+
        #         0.25*exp(logpdf(Normal(prop[k],0.0001),rho[k][ianim]))
        #         )
        #     else
        #         Loglike[k] -= log(prob_lim/3.0)
        #     end
        #
        #     if (rho[k][ianim]>0.0) && (rho[k][ianim]<1.0)
        #
        #         Loglike[k] += log(1.0-2.0*prob_lim/3.0)+ log(
        #         0.25*exp(logpdf(Normal(prop[k],0.1),rho[k][ianim]))+
        #         0.25*exp(logpdf(Normal(prop[k],0.01),rho[k][ianim]))+
        #         0.25*exp(logpdf(Normal(prop[k],0.001),rho[k][ianim]))+
        #         0.25*exp(logpdf(Normal(prop[k],0.0001),rho[k][ianim]))
        #         )
        #     else
        #         Loglike[k] += log(prob_lim/3.0)
        #     end
        #
        #     #prop[k] = rho[k][ianim]
        # end


        # for k in 1:kmax
        #     prop[k] = rho[k][ianim]+0.00001
        # end
        Dpsirho_prop     = deepcopy(Dpsirho)
        Drho_prop        = deepcopy(Drho)

        for k in 1:kmax
            WW = (ianim-1)*2 .+ [1;2]
            Dpsirho_prop[k][WW,WW]     = (1.0-prop[k]).*psi[k][ianim].*Matrix(I,2,2)
            Drho_prop[k][WW,WW]        = prop[k].*Matrix(I,2,2)
        end


        Cangle      = deepcopy(StartAngle[1])
        MatR        = zeros(Float64,nc,nc)
        MatRprop    = zeros(Float64,nc,nc)
        MatR2        = zeros(Float64,nc,nc)
        MatR2prop    = zeros(Float64,nc,nc)
        for i in 2:nt
            k       = zeta[i-1]
            yt1P    = Obs[i]  #view(Likelihood.data.data,2,:)
            yt      = Obs[i-1]

            for ianim2 in 1:nanim
                W = [1,2] .+(ianim2-1)*2
                MatR[W,W]       = [ cos(Cangle[ianim2]*rho[k][ianim2]) -sin(Cangle[ianim2]*rho[k][ianim2]); sin(Cangle[ianim2]*rho[k][ianim2])  cos(Cangle[ianim2]*rho[k][ianim2])  ]
                MatR2[W,W]       = [ cos(Cangle[ianim2] ) -sin(Cangle[ianim2] ); sin(Cangle[ianim2] )  cos(Cangle[ianim2] )  ]
                if ianim2!=ianim
                    MatR_prop[W,W]  = [ cos(Cangle[ianim2]*rho[k][ianim2]) -sin(Cangle[ianim2]*rho[k][ianim2]); sin(Cangle[ianim2]*rho[k][ianim2])  cos(Cangle[ianim2]*rho[k][ianim2])  ]
                    MatR2_prop[W,W]  = [ cos(Cangle[ianim2] ) -sin(Cangle[ianim2] ); sin(Cangle[ianim2] )  cos(Cangle[ianim2] )  ]
                else
                    MatR_prop[W,W]  = [ cos(Cangle[ianim2]*prop[k]) -sin(Cangle[ianim2]*prop[k]); sin(Cangle[ianim2]*prop[k])  cos(Cangle[ianim2]*prop[k])  ]
                    MatR2_prop[W,W]  = [ cos(Cangle[ianim2] ) -sin(Cangle[ianim2] ); sin(Cangle[ianim2] )  cos(Cangle[ianim2] )  ]
                end

            end

            ##
            Mean_prop   = transpose(MatR_prop)*(yt+Dpsirho_prop[k]*(mu0[k]-yt) + Drho_prop[k]*MatR2_prop*muC[k])
            Mean_acc    = transpose(MatR)*(yt+Dpsirho[k]*(mu0[k]-yt) + Drho[k]*MatR2*muC[k])


            Loglike[k] += logpdf(MvNormal(Mean_prop, sigma[k]), transpose(MatR_prop)*yt1P)
            Loglike[k] -= logpdf(MvNormal(Mean_acc,  sigma[k]),   transpose(MatR)*yt1P)
            ##

            for j in 1:nanim
                WW = [1 2] .+ (j-1)*2
                Cangle[j] = atan(yt1P[WW][2]-yt[WW][2],yt1P[WW][1]-yt[WW][1])
            end
            #print(Loglike)
            #print("\n\n")

        end

        for k in 1:kmax
            u = rand(rng,Uniform(0.0,1.0))
            # print(exp(Loglike[k]))
            # print("\n")
            # print(prop[k])
            # print(" ")
            # print(rhoPar.parameteracc[k][ianim])
            # print("\n")
            if u<exp(Loglike[k])
                rhoPar.parameteracc[k][ianim]  = prop[k]
                rhoPar.parameterprop[k][ianim] = prop[k]
            end
        end

        # print("LogLike")
        # print(Loglike)
        # print("\n")
    end
    # print("RHO")
    # print(rhoPar.parameteracc)
    # print("\n")


    return nothing

end
