function sample_missing!(rng::MersenneTwister,Likelihood::AbstractLikelihood,Missing::AbstractMissing)
    error(string("sample_missing not defined for Likelihood ", typeof(Likelihood) ) )
end

function sample_missing!(rng::MersenneTwister,Likelihood::Likelihood_OUmodel,Missing::MissingDoNotUpdate)
    return nothing
end

function sample_missing!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel},Missing::MissingDoUpdate)
    #nanim       = Likelihood[1].data.nanimals
    nanim       = size(Likelihood)[1]
    for ian in 1:nanim
        #print(ian)
        sample_missing!(rng,Likelihood[ian],Likelihood[ian].miss)
    end
end

#####
function sample_missing!(rng::MersenneTwister,Likelihood::Likelihood_OU_CircLinmodel,Missing::MissingDoUpdate)

    ### missing
    #Likelihood = MCMCLikelihood

    #print("dosamp")
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



    IndCol = Likelihood.miss.indexcol
    IndRow = Likelihood.miss.indexrow

    nmiss = size(IndRow,1)

    Cangle          = zeros(Float64,nanim)
    Cangle_prop          = zeros(Float64,nanim)
    MatR           = zeros(Float64,nc,nc)
    MatR_prop        = zeros(Float64,nc,nc)
    MatR2           = zeros(Float64,nc,nc)
    MatR2_prop        = zeros(Float64,nc,nc)
    for i in 1:nmiss

        LogLike     = 0.0

        irow        = IndRow[i]
        nmiss_i     = size(IndCol[i],1)
        IndexA      = IndCol[i]
        CondMean    = zeros(Float64,nmiss_i)

        CondVarInv  = zeros(Float64,nmiss_i,nmiss_i)

        Obs_prop    = deepcopy(Obs[irow])
        Obs_acc     = Obs[irow]

        if irow !=1#
            k    = Likelihood.clusterization.zeta[irow-1]

            if irow == 2
                Cangle      = deepcopy(StartAngle[1])
            else
                for j in 1:nanim
                    WW = [1 2] .+ (j-1)*2
                    Cangle[j] = atan(Obs[irow-1][WW[2]]-Obs[irow-2][WW[2]],Obs[irow-1][WW[1]]-Obs[irow-2][WW[1]])
                end
            end

            for ianim in 1:nanim
                W = [1,2] .+(ianim-1)*2
                MatR[W,W] = [ cos(Cangle[ianim]*rho[k][ianim]) -sin(Cangle[ianim]*rho[k][ianim]); sin(Cangle[ianim]*rho[k][ianim])  cos(Cangle[ianim]*rho[k][ianim])  ]
                MatR2[W,W] = [ cos(Cangle[ianim]) -sin(Cangle[ianim]); sin(Cangle[ianim])  cos(Cangle[ianim])  ]
            end





            yt         = Obs[irow-1]
            ytP1       = Obs[irow]

            # Mean_prop  = transpose(MatR)*(yt+Dpsirho[k]*(mu0[k]-yt) + Drho[k]*MatR*muC[k])
            #
            #
            # CondMean1, CondVar1, CondVarInv1 = compute_CondMeanAndVariance_MvNormal(IndexA,Mean_prop, sigmainv[k].mat,transpose(MatR)*ytP1)
            # Obs_prop[IndexA] = MatR[IndexA,IndexA]*rand(MvNormal(CondMean1, PDMat( Symmetric(CondVar1))  ))

            Mean_prop  = yt+Dpsirho[k]*(mu0[k]-yt) + Drho[k]*MatR2*muC[k]


            CondMean1, CondVar1, CondVarInv1 = compute_CondMeanAndVariance_MvNormal(IndexA,Mean_prop, sigmainv[k].mat,ytP1)
            Obs_prop[IndexA] = rand(rng,MvNormal(CondMean1, PDMat( Symmetric(MatR[IndexA,IndexA]*CondVar1*transpose(MatR[IndexA,IndexA])))  ))
        else

            for ia in IndexA
                Obs_prop[ia] = rand(rng,Normal(0.0,1.0))
            end
        end

        ####
        if irow <= nt-1

            k          = Likelihood.clusterization.zeta[irow]

            if irow >= 2
                for j in 1:nanim
                    WW = [1 2] .+ (j-1)*2
                    Cangle[j] = atan(Obs_acc[WW[2]]-Obs[irow-1][WW[2]],Obs_acc[WW[1]]-Obs[irow-1][WW[1]])
                end
                for j in 1:nanim
                    WW = [1 2] .+ (j-1)*2
                    Cangle_prop[j] = atan(Obs_prop[WW[2]]-Obs[irow-1][WW[2]],Obs_prop[WW[1]]-Obs[irow-1][WW[1]])
                end
                for ianim in 1:nanim
                    W = [1,2] .+(ianim-1)*2
                    MatR[W,W] = [ cos(Cangle[ianim]*rho[k][ianim]) -sin(Cangle[ianim]*rho[k][ianim]); sin(Cangle[ianim]*rho[k][ianim])  cos(Cangle[ianim]*rho[k][ianim])  ]
                    MatR2[W,W] = [ cos(Cangle[ianim] ) -sin(Cangle[ianim] ); sin(Cangle[ianim] )  cos(Cangle[ianim] )  ]
                end
                for ianim in 1:nanim
                    W = [1,2] .+(ianim-1)*2
                    MatR_prop[W,W] = [ cos(Cangle_prop[ianim]*rho[k][ianim]) -sin(Cangle_prop[ianim]*rho[k][ianim]); sin(Cangle_prop[ianim]*rho[k][ianim])  cos(Cangle_prop[ianim]*rho[k][ianim])  ]
                    MatR2_prop[W,W] = [ cos(Cangle_prop[ianim] ) -sin(Cangle_prop[ianim] ); sin(Cangle_prop[ianim] )  cos(Cangle_prop[ianim] )  ]
                end
            else

                Cangle              = deepcopy(StartAngle[1])
                Cangle_prop         = deepcopy(StartAngle[1])

                for ianim in 1:nanim
                    W = [1,2] .+(ianim-1)*2
                    MatR[W,W] = [ cos(Cangle[ianim]*rho[k][ianim]) -sin(Cangle[ianim]*rho[k][ianim]); sin(Cangle[ianim]*rho[k][ianim])  cos(Cangle[ianim]*rho[k][ianim])  ]
                    MatR2[W,W] = [ cos(Cangle[ianim]) -sin(Cangle[ianim] ); sin(Cangle[ianim])  cos(Cangle[ianim])  ]
                end
                for ianim in 1:nanim
                    W = [1,2] .+(ianim-1)*2
                    MatR_prop[W,W] = [ cos(Cangle_prop[ianim]*rho[k][ianim]) -sin(Cangle_prop[ianim]*rho[k][ianim]); sin(Cangle_prop[ianim]*rho[k][ianim])  cos(Cangle_prop[ianim]*rho[k][ianim])  ]
                    MatR2_prop[W,W] = [ cos(Cangle_prop[ianim] ) -sin(Cangle_prop[ianim] ); sin(Cangle_prop[ianim] )  cos(Cangle_prop[ianim] )  ]
                end
            end



            ytP1       = Obs[irow+1]
            Mean_acc   = transpose(MatR)*(Obs_acc+Dpsirho[k]*(mu0[k]-Obs_acc) + Drho[k]*MatR2*muC[k])
            Mean_prop  = transpose(MatR_prop)*(Obs_prop+Dpsirho[k]*(mu0[k]-Obs_prop) + Drho[k]*MatR2_prop*muC[k])

            LogLike   += logpdf(MvNormal(Mean_prop, sigma[k]), transpose(MatR_prop)*ytP1)
            LogLike   -= logpdf(MvNormal(Mean_acc, sigma[k]), transpose(MatR)*ytP1)
        end
        ####
        if irow <= nt-2

            k          = Likelihood.clusterization.zeta[irow+1]
            for j in 1:nanim
                WW = [1 2] .+ (j-1)*2
                Cangle[j] = atan(Obs[irow+1][WW[2]]-Obs_acc[WW[2]],Obs[irow+1][WW[1]]-Obs_acc[WW[1]])
            end
            for j in 1:nanim
                WW = [1 2] .+ (j-1)*2
                Cangle_prop[j] = atan(Obs[irow+1][WW[2]]-Obs_prop[WW[2]],Obs[irow+1][WW[1]]-Obs_prop[WW[1]])
            end
            for ianim in 1:nanim
                W = [1,2] .+(ianim-1)*2
                MatR[W,W] = [ cos(Cangle[ianim]*rho[k][ianim]) -sin(Cangle[ianim]*rho[k][ianim]); sin(Cangle[ianim]*rho[k][ianim])  cos(Cangle[ianim]*rho[k][ianim])  ]
                MatR2[W,W] = [ cos(Cangle[ianim] ) -sin(Cangle[ianim] ); sin(Cangle[ianim] )  cos(Cangle[ianim] )  ]
            end
            for ianim in 1:nanim
                W = [1,2] .+(ianim-1)*2
                MatR_prop[W,W] = [ cos(Cangle_prop[ianim]*rho[k][ianim]) -sin(Cangle_prop[ianim]*rho[k][ianim]); sin(Cangle_prop[ianim]*rho[k][ianim])  cos(Cangle_prop[ianim]*rho[k][ianim])  ]
                MatR2_prop[W,W] = [ cos(Cangle_prop[ianim] ) -sin(Cangle_prop[ianim] ); sin(Cangle_prop[ianim] )  cos(Cangle_prop[ianim] )  ]
            end

            ytP1       = Obs[irow+2]
            yt         = Obs[irow+1]

            Mean_acc   = transpose(MatR)*(yt+Dpsirho[k]*(mu0[k]-yt) + Drho[k]*MatR2*muC[k])
            Mean_prop  = transpose(MatR_prop)*(yt+Dpsirho[k]*(mu0[k]-yt) + Drho[k]*MatR2_prop*muC[k])

            LogLike   += logpdf(MvNormal(Mean_prop, sigma[k]), transpose(MatR_prop)*ytP1)
            LogLike   -= logpdf(MvNormal(Mean_acc, sigma[k]), transpose(MatR)*ytP1)
        end

        u = rand(rng,Uniform(0.0,1.0))

        if u<exp(LogLike)
            Obs_acc[IndexA] = Obs_prop[IndexA]
        end

    end

    #### angle
    for j in 1:nanim
        k    = Likelihood.clusterization.zeta[1]

        Angle_Acc   = deepcopy(StartAngle[1])
        Angle_prop  = deepcopy(StartAngle[1])

        Angle_prop[j] = rand(rng,Uniform(-pi, pi))


        LogLike = 0.0

        for ianim in 1:nanim
            W = [1,2] .+(ianim-1)*2
            MatR[W,W] = [ cos(Angle_Acc[ianim]*rho[k][ianim]) -sin(Angle_Acc[ianim]*rho[k][ianim]); sin(Angle_Acc[ianim]*rho[k][ianim])  cos(Angle_Acc[ianim]*rho[k][ianim])  ]
            MatR2[W,W] = [ cos(Angle_Acc[ianim] ) -sin(Angle_Acc[ianim] ); sin(Angle_Acc[ianim] )  cos(Angle_Acc[ianim] )  ]
        end

        for ianim in 1:nanim
            W = [1,2] .+(ianim-1)*2
            MatR_prop[W,W] = [ cos(Angle_prop[ianim]*rho[k][ianim]) -sin(Angle_prop[ianim]*rho[k][ianim]); sin(Angle_prop[ianim]*rho[k][ianim])  cos(Angle_prop[ianim]*rho[k][ianim])  ]
            MatR2_prop[W,W] = [ cos(Angle_prop[ianim] ) -sin(Angle_prop[ianim] ); sin(Angle_prop[ianim] )  cos(Angle_prop[ianim] )  ]
        end


        ytP1       = Obs[2]
        yt         = Obs[1]

        Mean_acc   = transpose(MatR)*(yt+Dpsirho[k]*(mu0[k]-yt) + Drho[k]*MatR2*muC[k])
        Mean_prop  = transpose(MatR_prop)*(yt+Dpsirho[k]*(mu0[k]-yt) + Drho[k]*MatR2_prop*muC[k])

        LogLike   += logpdf(MvNormal(Mean_prop, sigma[k]), transpose(MatR_prop)*ytP1)
        LogLike   -= logpdf(MvNormal(Mean_acc, sigma[k]), transpose(MatR)*ytP1)

        u = rand(rng,Uniform(0.0,1.0))

        if u<exp(LogLike)
            StartAngle[1][j] = Angle_prop[j]
        end
    end



    return nothing

end
