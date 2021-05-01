#### #### #### #### #### #### #### #### #### ####
#### #### #### sample mu0
#### #### #### #### #### #### #### #### #### ####

function sample_mu0!(rng::MersenneTwister,Likelihood::AbstractLikelihood, mu0Par::AbstractVecPar)
    error(string("sample_mu0 not defined for Likelihood ", typeof(Likelihood), " and mu0 ", typeof(mu0Par)) )
end





###

function sample_mu0!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, mu0Par::VecParMvNormal,parHier::OptionHierarchicalParameters)

    nanim_v2    = size(Likelihood)[1]
    kmax        = Likelihood[1].kmax

    InvMat_P      = Vector{Matrix{Float64}}()
    for k in 1:kmax
        push!(InvMat_P, deepcopy(mu0Par.prior_invmat.mat))
    end

    Mean_P      = Vector{Vector{Float64}}()
    for k in 1:kmax
        push!(Mean_P, deepcopy(mu0Par.prior.μ))
    end

    for ian in 1:nanim_v2

        nc          = Likelihood[ian].data.ncol
        nanim       = Likelihood[ian].data.nanimals
        nt          = Likelihood[ian].data.nt
        mu0         = Likelihood[ian].mu.parameteracc
        psi         = Likelihood[ian].nu.parameteracc
        rho         = Likelihood[ian].rho.parameteracc
        muC         = Likelihood[ian].eta.parameteracc
        sigmainv    = Likelihood[ian].sigma.parameteraccinv

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

        IdMatrix = Matrix(1.0*I, nc, nc)
        appVar =  Vector{Matrix{Float64}}()
        appMean =  Vector{Matrix{Float64}}()



        Cangle  = deepcopy(StartAngle[1])
        MatR = zeros(Float64,nc,nc)
        MatR2 = zeros(Float64,nc,nc)
        for i in 2:nt
            k       = zeta[i-1]
            yt1P    = Obs[i]  #view(Likelihood.data.data,2,:)
            yt      = Obs[i-1]

            for ianim in 1:nanim
                W = [1,2] .+(ianim-1)*2
                MatR[W,W] = [ cos(Cangle[ianim]*rho[k][ianim]) -sin(Cangle[ianim]*rho[k][ianim]); sin(Cangle[ianim]*rho[k][ianim])  cos(Cangle[ianim]*rho[k][ianim])  ]
                MatR2[W,W] = [ cos(Cangle[ianim] ) -sin(Cangle[ianim] ); sin(Cangle[ianim] )  cos(Cangle[ianim] )  ]
            end


            ##
            app         =  Dpsirho[k]*MatR*sigmainv[k].mat*transpose(MatR)
            InvMat_P[parHier.h_mu.clust[k,ian]] += app*Dpsirho[k]
            Mean_P[parHier.h_mu.clust[k,ian]]   += app*(yt1P-yt+Dpsirho[k]*yt-Drho[k]*MatR2*muC[k]  )
            ##

            for j in 1:nanim
                WW = [1 2] .+ (j-1)*2
                Cangle[j] = atan(yt1P[WW][2]-yt[WW][2],yt1P[WW][1]-yt[WW][1])
            end


        end
    end



    for k in 1:kmax
        Covmat = inv(PDMat(Symmetric(InvMat_P[k])))
        #PDMat(Symmetric(InvMat_P[k])).mat-InvMat_P[k]
        #B = PDMat(Symmetric(InvMat_P[k])).mat
        Mean   = Covmat*Mean_P[k]

        parHier.h_mu.par.parameteracc[k] = rand(rng,MvNormal(Mean,Covmat))
        parHier.h_mu.par.parameterprop[k] = deepcopy(mu0Par.parameteracc[k])
        #mu0Par.parameteracc[k] = rand(MvNormal(Mean,Covmat))
        #mu0Par.parameterprop[k] = deepcopy(mu0Par.parameteracc[k])
    end

    for ian in 1:nanim_v2
        mu0         = Likelihood[ian].mu
        for k in 1:kmax
            mu0.parameteracc[k]  = deepcopy(parHier.h_mu.par.parameteracc[parHier.h_mu.clust[k,ian]])
            mu0.parameterprop[k] = deepcopy(mu0.parameteracc[k])
        end
    end

    #isposdef(Likelihood.sigma.parameteraccinv[k])

    return nothing

end




###
function sample_mu0!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, mu0Par::VecParMvNormal)

    nanim = size(Likelihood)[1]
    for ian in 1:nanim
        sample_mu0!(rng,Likelihood[ian], Likelihood[ian].mu)
    end

end
function sample_mu0!(rng::MersenneTwister,Likelihood::Likelihood_OU_CircLinmodel, mu0Par::VecParMvNormal)

    # Likelihood  = MCMCLikelihood
    # mu0Par      = MCMCLikelihood.mu0
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    psi         = Likelihood.nu.parameteracc
    rho         = Likelihood.rho.parameteracc
    muC         = Likelihood.eta.parameteracc
    sigmainv    = Likelihood.sigma.parameteraccinv

    zeta        = Likelihood.clusterization.zeta
    Obs         = Likelihood.data.data

    StartAngle  = Likelihood.Angle.parameteracc

    InvMat_P      = Vector{Matrix{Float64}}()
    for k in 1:kmax
        push!(InvMat_P, deepcopy(mu0Par.prior_invmat.mat))
    end

    Mean_P      = Vector{Vector{Float64}}()
    for k in 1:kmax
        push!(Mean_P, deepcopy(mu0Par.prior.μ))
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

    IdMatrix = Matrix(1.0*I, nc, nc)
    appVar =  Vector{Matrix{Float64}}()
    appMean =  Vector{Matrix{Float64}}()



    Cangle  = deepcopy(StartAngle[1])
    MatR = zeros(Float64,nc,nc)
    MatR2 = zeros(Float64,nc,nc)
    for i in 2:nt
        k       = zeta[i-1]
        yt1P    = Obs[i]  #view(Likelihood.data.data,2,:)
        yt      = Obs[i-1]

        for ianim in 1:nanim
            W = [1,2] .+(ianim-1)*2
            MatR[W,W] = [ cos(Cangle[ianim]*rho[k][ianim]) -sin(Cangle[ianim]*rho[k][ianim]); sin(Cangle[ianim]*rho[k][ianim])  cos(Cangle[ianim]*rho[k][ianim])  ]
            MatR2[W,W] = [ cos(Cangle[ianim] ) -sin(Cangle[ianim] ); sin(Cangle[ianim] )  cos(Cangle[ianim] )  ]
        end

        ##
        app         =  Dpsirho[k]*MatR*sigmainv[k].mat*transpose(MatR)
        InvMat_P[k] += app*Dpsirho[k]
        Mean_P[k]   += app*(yt1P-yt+Dpsirho[k]*yt-Drho[k]*MatR2*muC[k]  )
        ##

        for j in 1:nanim
            WW = [1 2] .+ (j-1)*2
            Cangle[j] = atan(yt1P[WW][2]-yt[WW][2],yt1P[WW][1]-yt[WW][1])
        end


    end

    for k in 1:kmax
        Covmat = inv(PDMat(Symmetric(InvMat_P[k])))
        #PDMat(Symmetric(InvMat_P[k])).mat-InvMat_P[k]
        #B = PDMat(Symmetric(InvMat_P[k])).mat
        Mean   = Covmat*Mean_P[k]

        mu0Par.parameteracc[k] = rand(rng,MvNormal(Mean,Covmat))
        mu0Par.parameterprop[k] = deepcopy(mu0Par.parameteracc[k])
    end
    #isposdef(Likelihood.sigma.parameteraccinv[k])

    return nothing

end
