#
# function compute_logpdf(Likelihood::Likelihood_OUmodel, k::Int16, i::Int16)
#
#     nanim       = Likelihood.data.nanimals
# #log(sigmaPar.chol[1,1])
#
#     yP1         = Likelihood.data.data[i]
#     y           = Likelihood.data.data[i-1]
#     mu0Par      = Likelihood.mu0.parameteracc[k]
#     sigmaPar    = Likelihood.sigma.parameteracc[k]
#     psiPar      = Likelihood.psi.parameteracc[k]
#     #sigmaParInv    = Likelihood.sigma.parameteraccinv[k]
#
#     app_psi = repeat([1:nanim;], inner=2, outer=1)
#     Dpsi  =  diagm(psiPar[ app_psi ])
#
#     Mean = mu0Par+Dpsi*(y-mu0Par  )
#
#
#     # LogDet = Float64(0.0)
#     # for i in 1:sigmaPar.dim
#     #     LogDet  += log(sigmaPar.chol.L[i,i])
#     # end
#     # LogDet = 2.0*LogDet
#
#     #return -0.5*LogDet-0.5*transpose(yP1 .-Mean)*sigmaParInv.mat*(yP1 .-Mean)
#     return logpdf(MvNormal(Mean,sigmaPar),yP1)
# end
#


function compute_logpdf(Likelihood::Likelihood_OU_CircLinmodel, k::Int16, i::Int16)

    nanim       = Likelihood.data.nanimals
#log(sigmaPar.chol[1,1])
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt
    @inbounds yP1         = Likelihood.data.data[i]
    @inbounds y           = Likelihood.data.data[i-1]

    psi         = Likelihood.nu.parameteracc
    rho         = Likelihood.rho.parameteracc
    muC         = Likelihood.eta.parameteracc
    mu0         = Likelihood.mu.parameteracc
    sigma       = Likelihood.sigma.parameteracc

    zeta        = Likelihood.clusterization.zeta
    Obs         = Likelihood.data.data

    # print("Loke \n")
    # print(Likelihood.mu.parameteracc[k])
    #Dpsirho =  Matrix{Float64}()

    app         = repeat([1:nanim;], inner=2, outer=1)
    @inbounds Dpsirho   = diagm(  psi[k][ app ].*(1.0 .-rho[k][ app ])   )

    #Drho =  Matrix{Float64}()

    app  = repeat([1:nanim;], inner=2, outer=1)
    @inbounds Drho = diagm(rho[k][ app ]  )



    StartAngle  = Likelihood.Angle.parameteracc
    MatR = zeros(Float64,nc,nc)
    MatR2 = zeros(Float64,nc,nc)

    Cangle  = deepcopy(StartAngle[1])
    if i!=2
        for j in 1:nanim
            WW = [1 2] .+ (j-1)*2
            @inbounds Cangle[j] = atan(Obs[i-1][WW][2]-Obs[i-2][WW][2],Obs[i-1][WW][1]-Obs[i-2][WW][1])
        end
    end


    for ianim in 1:nanim
        W = [1,2] .+(ianim-1)*2
        @inbounds MatR[W,W] = [ cos(Cangle[ianim]*rho[k][ianim]) -sin(Cangle[ianim]*rho[k][ianim]); sin(Cangle[ianim]*rho[k][ianim])  cos(Cangle[ianim]*rho[k][ianim])  ]
        @inbounds MatR2[W,W] = [ cos(Cangle[ianim] ) -sin(Cangle[ianim] ); sin(Cangle[ianim] )  cos(Cangle[ianim] )  ]
    end

    @inbounds Mean = transpose(MatR)*(y+Dpsirho*(mu0[k]-y) + Drho*MatR2*muC[k])

    return logpdf(MvNormal(Mean,sigma[k]),transpose(MatR)*yP1)
end




function compute_logpdf(Likelihood_Joint::Vector{Likelihood_OU_CircLinmodel}, k::Int16, i::Int16)

    nanim_v2    = size(Likelihood_Joint)[1]
    logdens_joint   = Float64(0.0)
    for ian in 1:nanim_v2
        Likelihood  = Likelihood_Joint[ian]
        nanim       = Likelihood.data.nanimals
    #log(sigmaPar.chol[1,1])
        nc          = Likelihood.data.ncol
        nanim       = Likelihood.data.nanimals
        nt          = Likelihood.data.nt
        yP1         = Likelihood.data.data[i]
        y           = Likelihood.data.data[i-1]

        psi         = Likelihood.nu.parameteracc
        rho         = Likelihood.rho.parameteracc
        muC         = Likelihood.eta.parameteracc
        mu0         = Likelihood.mu.parameteracc
        sigma       = Likelihood.sigma.parameteracc

        zeta        = Likelihood.clusterization.zeta
        Obs         = Likelihood.data.data


        #Dpsirho =  Matrix{Float64}()

        app         = repeat([1:nanim;], inner=2, outer=1)
        Dpsirho   = diagm(  psi[k][ app ].*(1.0 .-rho[k][ app ])   )

        #Drho =  Matrix{Float64}()

        app  = repeat([1:nanim;], inner=2, outer=1)
        Drho = diagm(rho[k][ app ]  )



        StartAngle  = Likelihood.Angle.parameteracc
        MatR = zeros(Float64,nc,nc)
        MatR2 = zeros(Float64,nc,nc)

        Cangle  = deepcopy(StartAngle[1])
        if i!=2
            for j in 1:nanim
                WW = [1 2] .+ (j-1)*2
                Cangle[j] = atan(Obs[i-1][WW][2]-Obs[i-2][WW][2],Obs[i-1][WW][1]-Obs[i-2][WW][1])
            end
        end


        for ianim in 1:nanim
            W = [1,2] .+(ianim-1)*2
            MatR[W,W] = [ cos(Cangle[ianim]*rho[k][ianim]) -sin(Cangle[ianim]*rho[k][ianim]); sin(Cangle[ianim]*rho[k][ianim])  cos(Cangle[ianim]*rho[k][ianim])  ]
            MatR2[W,W] = [ cos(Cangle[ianim] ) -sin(Cangle[ianim] ); sin(Cangle[ianim] )  cos(Cangle[ianim] )  ]
        end

        Mean = transpose(MatR)*(y+Dpsirho*(mu0[k]-y) + Drho*MatR2*muC[k])

        logdens_joint += logpdf(MvNormal(Mean,sigma[k]),transpose(MatR)*yP1)
    end




    return logdens_joint
end
