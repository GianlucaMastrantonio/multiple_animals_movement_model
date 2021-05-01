#### #### #### #### #### #### #### #### #### ####
#### #### #### sample mu0
#### #### #### #### #### #### #### #### #### ####

function sample_psi!(rng::MersenneTwister,Likelihood::AbstractLikelihood, psiPar::AbstractVecPar)
    error(string("sample_psi not defined for Likelihood ", typeof(Likelihood), " and psi ", typeof(psiPar)) )
end

function sample_psi!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, psiPar::VecParMvUniform)

    nanim = size(Likelihood)[1]
    for ian in 1:nanim
        sample_psi!(rng,Likelihood[ian], Likelihood[ian].nu)
    end

end


function sample_psi!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, psiPar::VecParMvUniform, parHier::OptionHierarchicalParameters)

    # Likelihood  = MCMCLikelihood
    # psiPar      = MCMCLikelihood.psi

    nanim_v2    = size(Likelihood)[1]
    kmax        = Likelihood[1].kmax

    Mean_P = zeros(Float64,kmax)
    Var_P  = zeros(Float64,kmax)
    for ian in 1:nanim_v2

        nc          = Likelihood[ian].data.ncol
        nanim       = Likelihood[ian].data.nanimals
        nt          = Likelihood[ian].data.nt

        psi         = Likelihood[ian].nu.parameteracc
        rho         = Likelihood[ian].rho.parameteracc
        muC         = Likelihood[ian].eta.parameteracc
        mu0         = Likelihood[ian].mu.parameteracc
        sigmainv    = Likelihood[ian].sigma.parameteraccinv

        zeta        = Likelihood[ian].clusterization.zeta
        Obs         = Likelihood[ian].data.data

        StartAngle  = Likelihood[ian].Angle.parameteracc

        CovInv_AgivenB      = Vector{Matrix{Float64}}()
        for k in 1:kmax
            push!(CovInv_AgivenB, zeros(Float64,2,2))
        end

        CovABInvCovB      = Vector{Matrix{Float64}}()
        for k in 1:kmax
            push!(CovABInvCovB, zeros(Float64,2,nc-2))
        end

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
        MatR2 = zeros(Float64,nc,nc)

        for ianim in 1:nanim



            IndexA = Int32[((ianim-1)*2+1):((ianim-1)*2+2);]
            IndexB = Int32[1:nc;]
            IndexB = deleteat!(IndexB,IndexA)
            for k in 1:kmax
                CovABInvCovB[k], App ,CovInv_AgivenB[k]  = compute_AppCondMeanAndVariance_MvNormal(IndexA,Likelihood[ian].sigma.parameteraccinv[k].mat )
            end
            Cangle  = deepcopy(StartAngle[1])
            for i in 2:nt

                k       = zeta[i-1]
                yt1P    = Obs[i]  #view(Likelihood.data.data,2,:)
                yt      = Obs[i-1]

                for ianim in 1:nanim
                    W = [1,2] .+(ianim-1)*2
                    MatR[W,W] = [ cos(Cangle[ianim]*rho[k][ianim]) -sin(Cangle[ianim]*rho[k][ianim]); sin(Cangle[ianim]*rho[k][ianim])  cos(Cangle[ianim]*rho[k][ianim])  ]
                    MatR2[W,W] = [ cos(Cangle[ianim] ) -sin(Cangle[ianim] ); sin(Cangle[ianim] )  cos(Cangle[ianim] )  ]
                end

                Mean    = transpose(MatR)*(Obs[i-1]+Dpsirho[k]*(mu0[k]-Obs[i-1] )+Drho[k]*MatR2*muC[k])
                ObsP1R  = transpose(MatR)*Obs[i]
                ObsR    = transpose(MatR)*Obs[i-1]
                app     = ObsP1R-Mean
                PartialMeanCond = CovABInvCovB[k]*(app[IndexB])

                app1  = (ObsP1R-(Mean-transpose(MatR)*Dpsirho[k]*(mu0[k]-Obs[i-1] )))[IndexA] -PartialMeanCond

                app2 = (1.0-rho[k][ianim]) .* ( ( transpose(MatR)*mu0[k]- ObsR)[IndexA,:] )


                Mean_P[parHier.h_nu.clust[k,ian]]   += (transpose(app2)*CovInv_AgivenB[k]*app1)[1]
                Var_P[parHier.h_nu.clust[k,ian]]    += (transpose(app2)*CovInv_AgivenB[k]*app2)[1]

                for j in 1:nanim
                    WW = [1 2] .+ (j-1)*2
                    Cangle[j] = atan(yt1P[WW][2]-yt[WW][2],yt1P[WW][1]-yt[WW][1])
                end

            end



        end

    end


    Var_P   = 1.0 ./Var_P
    Mean_P  = Var_P .*Mean_P

    for k in 1:kmax
        ianim = 1
        a = psiPar.prior.v[1].a
        b = psiPar.prior.v[1].b

        if abs(Var_P[k])  != Inf
            # if  Mean_P[k]>a && Mean_P[k]<b
            #
            #     psiPar.parameteracc[k][ianim] = rand(Truncated(Normal(Mean_P[k],sqrt(Var_P[k])),a,b))
            #
            # else
            #
            #     u    = rand(Uniform(0.0,1.0))
            #     prop = 0.0
            #     if u<0.25
            #         prop = rand(Normal(psiPar.parameteracc[k][ianim]-a, 1.0))
            #     elseif u<0.5
            #         prop = rand(Normal(psiPar.parameteracc[k][ianim]-a, 0.1))
            #     elseif u<0.75
            #         prop = rand(Normal(psiPar.parameteracc[k][ianim]-a, 0.001))
            #     else
            #         prop = rand(Normal(psiPar.parameteracc[k][ianim]-a, 0.0001))
            #     end
            #     psiPar.parameterprop[k][ianim] = a+rem(prop,b-a ,RoundDown)
            #
            #     #psiPar.parameterprop[k][ianim] = rand(Uniform(a,b))
            #
            #     alpha =  logpdf(Normal(Mean_P[k],sqrt(Var_P[k])), psiPar.parameterprop[k][ianim])-logpdf(Normal(Mean_P[k],sqrt(Var_P[k])), psiPar.parameteracc[k][ianim])
            #
            #     if rand(Uniform(0.0,1.0))<exp(alpha)
            #         psiPar.parameteracc[k][ianim] = psiPar.parameterprop[k][ianim]
            #     end
            # end

            parHier.h_nu.par.parameteracc[k][ianim] = rand(rng,Truncated(Normal(Mean_P[k],sqrt(Var_P[k])),a,b))


        else
            parHier.h_nu.par.parameteracc[k][ianim] = rand(rng,Uniform(a,b))
        end
        parHier.h_nu.par.parameterprop[k][ianim] = parHier.h_nu.par.parameteracc[k][ianim]

    end

    for ian in 1:nanim_v2
        psi         = Likelihood[ian].nu
        for k in 1:kmax
            psi.parameteracc[k]  = deepcopy(parHier.h_nu.par.parameteracc[parHier.h_nu.clust[k,ian]])
            psi.parameterprop[k] = deepcopy(psi.parameteracc[k])
        end
    end

    return nothing

end



function sample_psi!(rng::MersenneTwister,Likelihood::Likelihood_OU_CircLinmodel, psiPar::VecParMvUniform)

    # Likelihood  = MCMCLikelihood
    # psiPar      = MCMCLikelihood.psi

    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    psi         = Likelihood.nu.parameteracc
    rho         = Likelihood.rho.parameteracc
    muC         = Likelihood.eta.parameteracc
    mu0         = Likelihood.mu.parameteracc
    sigmainv    = Likelihood.sigma.parameteraccinv

    zeta        = Likelihood.clusterization.zeta
    Obs         = Likelihood.data.data

    StartAngle  = Likelihood.Angle.parameteracc

    CovInv_AgivenB      = Vector{Matrix{Float64}}()
    for k in 1:kmax
        push!(CovInv_AgivenB, zeros(Float64,2,2))
    end

    CovABInvCovB      = Vector{Matrix{Float64}}()
    for k in 1:kmax
        push!(CovABInvCovB, zeros(Float64,2,nc-2))
    end

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
    MatR2 = zeros(Float64,nc,nc)

    for ianim in 1:nanim

        Mean_P = zeros(Float64,kmax)
        Var_P  = zeros(Float64,kmax)

        IndexA = Int32[((ianim-1)*2+1):((ianim-1)*2+2);]
        IndexB = Int32[1:nc;]
        IndexB = deleteat!(IndexB,IndexA)
        for k in 1:kmax
            CovABInvCovB[k], App ,CovInv_AgivenB[k]  = compute_AppCondMeanAndVariance_MvNormal(IndexA,Likelihood.sigma.parameteraccinv[k].mat )
        end
        Cangle  = deepcopy(StartAngle[1])
        for i in 2:nt

            k       = zeta[i-1]
            yt1P    = Obs[i]  #view(Likelihood.data.data,2,:)
            yt      = Obs[i-1]

            for ianim in 1:nanim
                W = [1,2] .+(ianim-1)*2
                MatR[W,W] = [ cos(Cangle[ianim]*rho[k][ianim]) -sin(Cangle[ianim]*rho[k][ianim]); sin(Cangle[ianim]*rho[k][ianim])  cos(Cangle[ianim]*rho[k][ianim])  ]
                MatR2[W,W] = [ cos(Cangle[ianim] ) -sin(Cangle[ianim] ); sin(Cangle[ianim] )  cos(Cangle[ianim] )  ]
            end

            Mean    = transpose(MatR)*(Obs[i-1]+Dpsirho[k]*(mu0[k]-Obs[i-1] )+Drho[k]*MatR2*muC[k])
            ObsP1R  = transpose(MatR)*Obs[i]
            ObsR    = transpose(MatR)*Obs[i-1]
            app     = ObsP1R-Mean
            PartialMeanCond = CovABInvCovB[k]*(app[IndexB])

            app1  = (ObsP1R-(Mean-transpose(MatR)*Dpsirho[k]*(mu0[k]-Obs[i-1] )))[IndexA] -PartialMeanCond

            app2 = (1.0-rho[k][ianim]) .* ( ( transpose(MatR)*mu0[k]- ObsR)[IndexA,:] )


            Mean_P[k]   += (transpose(app2)*CovInv_AgivenB[k]*app1)[1]
            Var_P[k]    += (transpose(app2)*CovInv_AgivenB[k]*app2)[1]

            for j in 1:nanim
                WW = [1 2] .+ (j-1)*2
                Cangle[j] = atan(yt1P[WW][2]-yt[WW][2],yt1P[WW][1]-yt[WW][1])
            end

        end
        Var_P   = 1.0 ./Var_P
        Mean_P  = Var_P .*Mean_P

        for k in 1:kmax

            a = psiPar.prior.v[ianim].a
            b = psiPar.prior.v[ianim].b

            if abs(Var_P[k])  != Inf
                # if  Mean_P[k]>a && Mean_P[k]<b
                #
                #     psiPar.parameteracc[k][ianim] = rand(Truncated(Normal(Mean_P[k],sqrt(Var_P[k])),a,b))
                #
                # else
                #
                #     u    = rand(Uniform(0.0,1.0))
                #     prop = 0.0
                #     if u<0.25
                #         prop = rand(Normal(psiPar.parameteracc[k][ianim]-a, 1.0))
                #     elseif u<0.5
                #         prop = rand(Normal(psiPar.parameteracc[k][ianim]-a, 0.1))
                #     elseif u<0.75
                #         prop = rand(Normal(psiPar.parameteracc[k][ianim]-a, 0.001))
                #     else
                #         prop = rand(Normal(psiPar.parameteracc[k][ianim]-a, 0.0001))
                #     end
                #     psiPar.parameterprop[k][ianim] = a+rem(prop,b-a ,RoundDown)
                #
                #     #psiPar.parameterprop[k][ianim] = rand(Uniform(a,b))
                #
                #     alpha =  logpdf(Normal(Mean_P[k],sqrt(Var_P[k])), psiPar.parameterprop[k][ianim])-logpdf(Normal(Mean_P[k],sqrt(Var_P[k])), psiPar.parameteracc[k][ianim])
                #
                #     if rand(Uniform(0.0,1.0))<exp(alpha)
                #         psiPar.parameteracc[k][ianim] = psiPar.parameterprop[k][ianim]
                #     end
                # end

                psiPar.parameteracc[k][ianim] = rand(rng,Truncated(Normal(Mean_P[k],sqrt(Var_P[k])),a,b))


            else
                psiPar.parameteracc[k][ianim] = rand(rng,Uniform(a,b))
            end
            psiPar.parameterprop[k][ianim] = psiPar.parameteracc[k][ianim]

        end
        for k in 1:kmax
            for ip in 1:2
                Dpsirho[k][IndexA[ip],IndexA[ip]] = psiPar.parameteracc[k][ianim]*(1.0 .-rho[k][ianim])
            end
        end


    end



    return nothing

end
