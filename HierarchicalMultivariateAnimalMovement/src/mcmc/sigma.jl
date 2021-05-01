function sample_sigma!(rng::MersenneTwister,Likelihood::AbstractLikelihood, sigmaPar::AbstractPosDefMatPar)
    error(string("sample_sigma not defined for Likelihood ", typeof(Likelihood), " and sigma ", typeof(sigmaPar)) )
end

function sample_sigma!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, sigmaPar::PosDefMatInverseWishart)

    nanim = size(Likelihood)[1]
    for ian in 1:nanim
        sample_sigma!(rng,Likelihood[ian], Likelihood[ian].sigma)
    end

end


function sample_sigma!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, sigmaPar::PosDefMatInverseWishart,parHier::OptionHierarchicalParameters)

    #Likelihood  = MCMCLikelihood
    #sigmaPar      = MCMCLikelihood[1].sigma
    #parHier = HierarchicalParameters
    nanim_v2    = size(Likelihood)[1]
    kmax        = Likelihood[1].kmax

    Mat_P =  Vector{Matrix{Float64}}()
    for k in 1:kmax
        push!(Mat_P, deepcopy(sigmaPar.prior.Ψ.mat))
    end

    ParNum = ones(Float64,kmax)*sigmaPar.prior.df

    for ian in 1:nanim_v2

        nc          = Likelihood[ian].data.ncol
        nanim       = Likelihood[ian].data.nanimals
        nt          = Likelihood[ian].data.nt

        psi         = Likelihood[ian].nu.parameteracc
        rho         = Likelihood[ian].rho.parameteracc
        muC         = Likelihood[ian].eta.parameteracc
        mu0         = Likelihood[ian].mu.parameteracc
        #sigmainv    = Likelihood.sigma.parameteraccinv

        zeta        = Likelihood[ian].clusterization.zeta
        Obs         = Likelihood[ian].data.data

        StartAngle  = Likelihood[ian].Angle.parameteracc


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




        Cangle  = deepcopy(StartAngle[1])
        MatR    = zeros(Float64,nc,nc)
        MatR2    = zeros(Float64,nc,nc)
        for i in 2:nt
            k       = zeta[i-1]
            yt1P    = Obs[i]  #view(Likelihood.data.data,2,:)
            yt      = Obs[i-1]

            for ianim in 1:nanim
                W = [1,2] .+(ianim-1)*2
                MatR[W,W] = [ cos(Cangle[ianim]*rho[k][ianim]) -sin(Cangle[ianim]*rho[k][ianim]); sin(Cangle[ianim]*rho[k][ianim])  cos(Cangle[ianim]*rho[k][ianim])  ]
                MatR2[W,W] = [ cos(Cangle[ianim] ) -sin(Cangle[ianim] ); sin(Cangle[ianim] )  cos(Cangle[ianim] )  ]
            end

            SumsObs   = yt1P-(yt+Dpsirho[k]*(mu0[k]-yt) +Drho[k]*MatR2*muC[k])

            Mat_P[parHier.h_sigma.clust[k,ian]] += transpose(MatR)*SumsObs*transpose(SumsObs)*MatR
            ParNum[parHier.h_sigma.clust[k,ian]] += Float64(1)

            # if parHier.h_sigma.clust[k,ian] == 114
            #     app = Mat_P[parHier.h_sigma.clust[k,ian]]
            #     # print("Start \n",app[1,2]/(app[1,1]*app[2,2])^0.5)
            #     # print("\n")
            #     # print(app, "\n")
            #     # print(SumsObs,"\n")
            #     # print(MatR)
            # end
            for j in 1:nanim
                WW = [1 2] .+ (j-1)*2
                Cangle[j] = atan(yt1P[WW][2]-yt[WW][2],yt1P[WW][1]-yt[WW][1])
            end

        end
    end



    for k in 1:kmax
        #print(k)
        par1 = ParNum[k]
        parHier.h_sigma.par.parameteracc[k]    = PDMat(rand(rng,InverseWishart(par1,PDMat(Symmetric(Mat_P[k]))  )))
        #parHier.h_sigma.par.parameteraccinv[k] = inv(sigmaPar.parameteracc[k])
        parHier.h_sigma.par.parameteraccinv[k] = inv(parHier.h_sigma.par.parameteracc[k])

        parHier.h_sigma.par.parameterprop[k]    = deepcopy(parHier.h_sigma.par.parameteracc[k])
        parHier.h_sigma.par.parameterpropinv[k] = deepcopy(parHier.h_sigma.par.parameteraccinv[k])
    end
    for ian in 1:nanim_v2
        sigma         = Likelihood[ian].sigma
        for k in 1:kmax
            sigma.parameteracc[k]  = deepcopy(parHier.h_sigma.par.parameteracc[parHier.h_sigma.clust[k,ian]])
            sigma.parameterprop[k] = deepcopy(sigma.parameteracc[k])

            sigma.parameteraccinv[k]  = deepcopy(parHier.h_sigma.par.parameteraccinv[parHier.h_sigma.clust[k,ian]])
            sigma.parameterpropinv[k] = deepcopy(sigma.parameteraccinv[k])
        end
    end
    return nothing

end



function sample_sigma!(rng::MersenneTwister,Likelihood::Likelihood_OU_CircLinmodel, sigmaPar::PosDefMatInverseWishart)

    # Likelihood  = MCMCLikelihood
    # sigmaPar      = MCMCLikelihood.sigma

    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    psi         = Likelihood.nu.parameteracc
    rho         = Likelihood.rho.parameteracc
    muC         = Likelihood.eta.parameteracc
    mu0         = Likelihood.mu.parameteracc
    #sigmainv    = Likelihood.sigma.parameteraccinv

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


    Mat_P =  Vector{Matrix{Float64}}()
    for k in 1:kmax
        push!(Mat_P, deepcopy(sigmaPar.prior.Ψ.mat))
    end

    Cangle  = deepcopy(StartAngle[1])
    MatR    = zeros(Float64,nc,nc)
    MatR2    = zeros(Float64,nc,nc)
    for i in 2:nt
        k       = zeta[i-1]
        yt1P    = Obs[i]  #view(Likelihood.data.data,2,:)
        yt      = Obs[i-1]

        for ianim in 1:nanim
            W = [1,2] .+(ianim-1)*2
            MatR[W,W] = [ cos(Cangle[ianim]*rho[k][ianim]) -sin(Cangle[ianim]*rho[k][ianim]); sin(Cangle[ianim]*rho[k][ianim])  cos(Cangle[ianim]*rho[k][ianim])  ]
            MatR2[W,W] = [ cos(Cangle[ianim] ) -sin(Cangle[ianim] ); sin(Cangle[ianim] )  cos(Cangle[ianim] )  ]
        end

        SumsObs   = yt1P-(yt+Dpsirho[k]*(mu0[k]-yt) +Drho[k]*MatR2*muC[k])

        Mat_P[k] += transpose(MatR)*SumsObs*transpose(SumsObs)*MatR

        for j in 1:nanim
            WW = [1 2] .+ (j-1)*2
            Cangle[j] = atan(yt1P[WW][2]-yt[WW][2],yt1P[WW][1]-yt[WW][1])
        end

    end


    for k in 1:kmax
        par1 = sigmaPar.prior.df+typeof(sigmaPar.prior.df)(Likelihood.clusterization.n_obsC[k])
        sigmaPar.parameteracc[k]    = PDMat(rand(rng,InverseWishart(par1,PDMat(Symmetric(Mat_P[k]))  )))
        sigmaPar.parameteraccinv[k] = inv(sigmaPar.parameteracc[k])

        sigmaPar.parameterprop[k] = deepcopy(sigmaPar.parameteracc[k])
        sigmaPar.parameterpropinv[k] = deepcopy(sigmaPar.parameteraccinv[k])
    end

    return nothing

end
