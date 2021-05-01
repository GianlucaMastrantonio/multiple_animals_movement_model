
struct MCMCutils{T<:Integer,Td<:AbstractDict}
    iterations::T
    burnin::T
    thin::T
    nsamplesave::T
    mcmc_out::Td
    indexsave::Vector{T}

    function MCMCutils{T,Td}(iterations::T,burnin::T,thin::T,nsamplesave::T,mcmc_out::Td) where {T<:Integer,Td<:AbstractDict}
        if iterations <= burnin error("iterations must be greater than burnin") end
        if iterations < 1 error("iterations must greater than zero") end
        if burnin < 1 error("burnin must greater than zero") end
        if thin < 1 error("thin must greater than zero")  end
        indexsave = ones(T(1))
        new{T,Td}(iterations,burnin,thin,nsamplesave,mcmc_out,indexsave)
    end

end

# function MCMCutils(iterations::T,burnin::T,thin::T) where {T<:Integer}
#
#     MCMCutils{T}(iterations,burnin,thin,T((iterations-burnin)/thin))
# end



##### HMM DIVISO - gerarchico
function OptionsMCMC(Likelihood::Likelihood_OU_CircLinmodel, Clusterization::Clusterization_HDPHMM, HierarchicalParameters::OptionHierarchicalParameters;
    iterations::T = 100,
    burnin::T = 50,
    thin::T = 10,

    ) where {T<:Integer}

    if iterations <= burnin error("iterations must be greater than burnin") end
    if iterations < 1 error("iterations must greater than zero") end
    if burnin < 1 error("burnin must greater than zero") end
    if thin < 1 error("thin must greater than zero")  end

    nsamplesave = T((iterations-burnin)/thin)
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    if Likelihood.SaveMissing == false
        mcmc_out = Dict{String,Union{Array{Int16,2},Array{Float64,3},Array{Int16,2}} }(
        "mu"   => Array{Float64,3}(undef,nsamplesave,nc,kmax),
        "nu"   => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
        "sigma" => Array{Float64,3}(undef,nsamplesave,nc*nc,kmax),
        "zeta"  => Array{Int16,2}(undef,nsamplesave,nt-1),
        "pi"    => Array{Float64,3}(undef,nsamplesave,kmax,kmax),
        "eta"   => Array{Float64,3}(undef,nsamplesave,nc,kmax),
        "rho"   => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
        "z_mu"   => Array{Int16,2}(undef,nsamplesave,kmax),
        "z_nu"   => Array{Int16,2}(undef,nsamplesave,kmax),
        "z_sigma" => Array{Int16,2}(undef,nsamplesave,kmax),
        "z_eta"   => Array{Int16,2}(undef,nsamplesave,kmax),
        "z_rho"   => Array{Int16,2}(undef,nsamplesave,kmax)
        )
    else
        mcmc_out = Dict{String,Union{Array{Int16,2},Array{Float64,3},Array{Int16,2},Array{Int32,1},  Array{Array{Int32,1},1}} }(
        "mu"   => Array{Float64,3}(undef,nsamplesave,nc,kmax),
        "nu"   => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
        "sigma" => Array{Float64,3}(undef,nsamplesave,nc*nc,kmax),
        "zeta"  => Array{Int16,2}(undef,nsamplesave,nt-1),
        "pi"    => Array{Float64,3}(undef,nsamplesave,kmax,kmax),
        "eta"   => Array{Float64,3}(undef,nsamplesave,nc,kmax),
        "rho"   => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
        "missing" => Array{Float64,3}(undef,nsamplesave,size(Likelihood.miss.indexrow,1),Likelihood.data.ncol),
        "missingIndexRow" => Likelihood.miss.indexrow,
        "missingIndexCol" => Likelihood.miss.indexcol,
        "z_mu"   => Array{Int16,2}(undef,nsamplesave,kmax),
        "z_nu"   => Array{Int16,2}(undef,nsamplesave,kmax),
        "z_sigma" => Array{Int16,2}(undef,nsamplesave,kmax),
        "z_eta"   => Array{Int16,2}(undef,nsamplesave,kmax),
        "z_rho"   => Array{Int16,2}(undef,nsamplesave,kmax)
        )


    end
 #Likelihood = MCMCLikelihood



    return MCMCutils{T,typeof(mcmc_out)}(iterations,burnin,thin,nsamplesave,mcmc_out)
end

function save_posteriorsamples!(MCMCout::Vector{MCMCutils},Likelihood::Vector{Likelihood_OU_CircLinmodel}, Clusterization::Vector{Clusterization_HDPHMM},HP::OptionHierarchicalParameters)

    nanim = size(Likelihood)[1]

    for ian in 1:nanim
        save_posteriorsamples!(MCMCout[ian],Likelihood[ian], Clusterization[ian],HP,Int16(ian))
    end

end

function save_posteriorsamples!(MCMCout::MCMCutils,Likelihood::Likelihood_OU_CircLinmodel, Clusterization::Clusterization_HDPHMM,HP::OptionHierarchicalParameters, nanim_v2::Int16)

    # Likelihood = MCMCLikelihood
    # Clusterization = MCMCClusterization

    indOut      = MCMCout.indexsave

    if indOut[1]>MCMCout.nsamplesave error("indOut>nsamplesave") end

    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    if Likelihood.SaveMissing == false
        for k in 1:kmax
            MCMCout.mcmc_out["mu"][indOut[1],:,k] = Likelihood.mu.parameteracc[k]
            MCMCout.mcmc_out["nu"][indOut[1],:,k] = Likelihood.nu.parameteracc[k]
            MCMCout.mcmc_out["sigma"][indOut[1],:,k] = vec(Likelihood.sigma.parameteracc[k].mat)
            MCMCout.mcmc_out["pi"][indOut[1],:,k] = Clusterization.pi.parameteracc[k]
            MCMCout.mcmc_out["zeta"][indOut[1],:] = Clusterization.clusterization.zeta[1:(nt-1)]
            MCMCout.mcmc_out["eta"][indOut[1],:,k] = Likelihood.eta.parameteracc[k]
            MCMCout.mcmc_out["rho"][indOut[1],:,k] = Likelihood.rho.parameteracc[k]


            MCMCout.mcmc_out["z_mu"][indOut[1],k] = HP.h_mu.clust[k,nanim_v2]
            MCMCout.mcmc_out["z_nu"][indOut[1],k] = HP.h_nu.clust[k,nanim_v2]
            MCMCout.mcmc_out["z_sigma"][indOut[1],k] = HP.h_sigma.clust[k,nanim_v2]
            MCMCout.mcmc_out["z_eta"][indOut[1],k] = HP.h_eta.clust[k,nanim_v2]
            MCMCout.mcmc_out["z_rho"][indOut[1],k] = HP.h_rho.clust[k,nanim_v2]
        end
    else
        for k in 1:kmax
            MCMCout.mcmc_out["mu"][indOut[1],:,k] = Likelihood.mu.parameteracc[k]
            MCMCout.mcmc_out["nu"][indOut[1],:,k] = Likelihood.nu.parameteracc[k]
            MCMCout.mcmc_out["sigma"][indOut[1],:,k] = vec(Likelihood.sigma.parameteracc[k].mat)
            MCMCout.mcmc_out["pi"][indOut[1],:,k] = Clusterization.pi.parameteracc[k]
            MCMCout.mcmc_out["zeta"][indOut[1],:] = Clusterization.clusterization.zeta[1:(nt-1)]
            MCMCout.mcmc_out["eta"][indOut[1],:,k] = Likelihood.eta.parameteracc[k]
            MCMCout.mcmc_out["rho"][indOut[1],:,k] = Likelihood.rho.parameteracc[k]

            MCMCout.mcmc_out["missing"][indOut[1],:,:]  =  Matrix(transpose(hcat(Likelihood.data.data[Likelihood.miss.indexrow,:]...)))

            MCMCout.mcmc_out["z_mu"][indOut[1],k] = HP.h_mu.clust[k,nanim_v2]
            MCMCout.mcmc_out["z_nu"][indOut[1],k] = HP.h_nu.clust[k,nanim_v2]
            MCMCout.mcmc_out["z_sigma"][indOut[1],k] = HP.h_sigma.clust[k,nanim_v2]
            MCMCout.mcmc_out["z_eta"][indOut[1],k] = HP.h_eta.clust[k,nanim_v2]
            MCMCout.mcmc_out["z_rho"][indOut[1],k] = HP.h_rho.clust[k,nanim_v2]

            # for irow in 1:size(Likelihood.miss.indexrow,1)
            #
            # end
        end
    end

    MCMCout.indexsave[1] += one(MCMCout.indexsave[1])
#Likelihood.data.data[Likelihood.miss.indexrow,:]

#Matrix(Likelihood.data.data[Likelihood.miss.indexrow,:])
    return nothing
end



#####
##### HMM DIVISO
function OptionsMCMC(Likelihood::Likelihood_OU_CircLinmodel, Clusterization::Clusterization_HDPHMM;
    iterations::T = 100,
    burnin::T = 50,
    thin::T = 10,

    ) where {T<:Integer}

    if iterations <= burnin error("iterations must be greater than burnin") end
    if iterations < 1 error("iterations must greater than zero") end
    if burnin < 1 error("burnin must greater than zero") end
    if thin < 1 error("thin must greater than zero")  end

    nsamplesave = T((iterations-burnin)/thin)
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    if Likelihood.SaveMissing == false
        mcmc_out = Dict{String,Union{Array{Float64,3},Array{Int16,2}} }(
        "mu"   => Array{Float64,3}(undef,nsamplesave,nc,kmax),
        "nu"   => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
        "sigma" => Array{Float64,3}(undef,nsamplesave,nc*nc,kmax),
        "zeta"  => Array{Int16,2}(undef,nsamplesave,nt-1),
        "pi"    => Array{Float64,3}(undef,nsamplesave,kmax,kmax),
        "eta"   => Array{Float64,3}(undef,nsamplesave,nc,kmax),
        "rho"   => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
        )
    else
        mcmc_out = Dict{String,Union{Array{Float64,3},Array{Int16,2},Array{Int32,1},  Array{Array{Int32,1},1}} }(
        "mu"   => Array{Float64,3}(undef,nsamplesave,nc,kmax),
        "nu"   => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
        "sigma" => Array{Float64,3}(undef,nsamplesave,nc*nc,kmax),
        "zeta"  => Array{Int16,2}(undef,nsamplesave,nt-1),
        "pi"    => Array{Float64,3}(undef,nsamplesave,kmax,kmax),
        "eta"   => Array{Float64,3}(undef,nsamplesave,nc,kmax),
        "rho"   => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
        "missing" => Array{Float64,3}(undef,nsamplesave,size(Likelihood.miss.indexrow,1),Likelihood.data.ncol),
        "missingIndexRow" => Likelihood.miss.indexrow,
        "missingIndexCol" => Likelihood.miss.indexcol,
        )


    end
 #Likelihood = MCMCLikelihood



    return MCMCutils{T,typeof(mcmc_out)}(iterations,burnin,thin,nsamplesave,mcmc_out)
end

function save_posteriorsamples!(MCMCout::Vector{MCMCutils},Likelihood::Vector{Likelihood_OU_CircLinmodel}, Clusterization::Vector{Clusterization_HDPHMM})

    nanim = size(Likelihood)[1]

    for ian in 1:nanim
        save_posteriorsamples!(MCMCout[ian],Likelihood[ian], Clusterization[ian])
    end

end

function save_posteriorsamples!(MCMCout::MCMCutils,Likelihood::Likelihood_OU_CircLinmodel, Clusterization::Clusterization_HDPHMM)

    # Likelihood = MCMCLikelihood
    # Clusterization = MCMCClusterization

    indOut      = MCMCout.indexsave

    if indOut[1]>MCMCout.nsamplesave error("indOut>nsamplesave") end

    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    if Likelihood.SaveMissing == false
        for k in 1:kmax
            MCMCout.mcmc_out["mu"][indOut[1],:,k] = Likelihood.mu.parameteracc[k]
            MCMCout.mcmc_out["nu"][indOut[1],:,k] = Likelihood.nu.parameteracc[k]
            MCMCout.mcmc_out["sigma"][indOut[1],:,k] = vec(Likelihood.sigma.parameteracc[k].mat)
            MCMCout.mcmc_out["pi"][indOut[1],:,k] = Clusterization.pi.parameteracc[k]
            MCMCout.mcmc_out["zeta"][indOut[1],:] = Clusterization.clusterization.zeta[1:(nt-1)]
            MCMCout.mcmc_out["eta"][indOut[1],:,k] = Likelihood.eta.parameteracc[k]
            MCMCout.mcmc_out["rho"][indOut[1],:,k] = Likelihood.rho.parameteracc[k]
        end
    else
        for k in 1:kmax
            MCMCout.mcmc_out["mu"][indOut[1],:,k] = Likelihood.mu.parameteracc[k]
            MCMCout.mcmc_out["nu"][indOut[1],:,k] = Likelihood.nu.parameteracc[k]
            MCMCout.mcmc_out["sigma"][indOut[1],:,k] = vec(Likelihood.sigma.parameteracc[k].mat)
            MCMCout.mcmc_out["pi"][indOut[1],:,k] = Clusterization.pi.parameteracc[k]
            MCMCout.mcmc_out["zeta"][indOut[1],:] = Clusterization.clusterization.zeta[1:(nt-1)]
            MCMCout.mcmc_out["eta"][indOut[1],:,k] = Likelihood.eta.parameteracc[k]
            MCMCout.mcmc_out["rho"][indOut[1],:,k] = Likelihood.rho.parameteracc[k]

            MCMCout.mcmc_out["missing"][indOut[1],:,:]  =  Matrix(transpose(hcat(Likelihood.data.data[Likelihood.miss.indexrow,:]...)))
            # for irow in 1:size(Likelihood.miss.indexrow,1)
            #
            # end
        end
    end

    MCMCout.indexsave[1] += one(MCMCout.indexsave[1])
#Likelihood.data.data[Likelihood.miss.indexrow,:]

#Matrix(Likelihood.data.data[Likelihood.miss.indexrow,:])
    return nothing
end





##### HMM JOINT
function OptionsMCMC(Likelihood::Likelihood_OU_CircLinmodel, Clusterization::Clusterization_HDPHMM_Joint;
    iterations::T = 100,
    burnin::T = 50,
    thin::T = 10,

    ) where {T<:Integer}

    if iterations <= burnin error("iterations must be greater than burnin") end
    if iterations < 1 error("iterations must greater than zero") end
    if burnin < 1 error("burnin must greater than zero") end
    if thin < 1 error("thin must greater than zero")  end

    nsamplesave = T((iterations-burnin)/thin)
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    if Likelihood.SaveMissing == false
        mcmc_out = Dict{String,Union{Array{Float64,3},Array{Int16,2}} }(
        "mu"   => Array{Float64,3}(undef,nsamplesave,nc,kmax),
        "nu"   => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
        "sigma" => Array{Float64,3}(undef,nsamplesave,nc*nc,kmax),
        "zeta"  => Array{Int16,2}(undef,nsamplesave,nt-1),
        "pi"    => Array{Float64,3}(undef,nsamplesave,kmax,kmax),
        "eta"   => Array{Float64,3}(undef,nsamplesave,nc,kmax),
        "rho"   => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
        )
    else
        mcmc_out = Dict{String,Union{Array{Float64,3},Array{Int16,2},Array{Int32,1},  Array{Array{Int32,1},1}} }(
        "mu"   => Array{Float64,3}(undef,nsamplesave,nc,kmax),
        "nu"   => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
        "sigma" => Array{Float64,3}(undef,nsamplesave,nc*nc,kmax),
        "zeta"  => Array{Int16,2}(undef,nsamplesave,nt-1),
        "pi"    => Array{Float64,3}(undef,nsamplesave,kmax,kmax),
        "eta"   => Array{Float64,3}(undef,nsamplesave,nc,kmax),
        "rho"   => Array{Float64,3}(undef,nsamplesave,nanim,kmax),
        "missing" => Array{Float64,3}(undef,nsamplesave,size(Likelihood.miss.indexrow,1),Likelihood.data.ncol),
        "missingIndexRow" => Likelihood.miss.indexrow,
        "missingIndexCol" => Likelihood.miss.indexcol,
        )


    end
 #Likelihood = MCMCLikelihood



    return MCMCutils{T,typeof(mcmc_out)}(iterations,burnin,thin,nsamplesave,mcmc_out)
end

function save_posteriorsamples!(MCMCout::Vector{MCMCutils},Likelihood::Vector{Likelihood_OU_CircLinmodel}, Clusterization::Vector{Clusterization_HDPHMM_Joint})

    nanim = size(Likelihood)[1]

    for ian in 1:nanim
        save_posteriorsamples!(MCMCout[ian],Likelihood[ian], Clusterization[ian])
    end

end

function save_posteriorsamples!(MCMCout::MCMCutils,Likelihood::Likelihood_OU_CircLinmodel, Clusterization::Clusterization_HDPHMM_Joint)

    # Likelihood = MCMCLikelihood
    # Clusterization = MCMCClusterization

    indOut      = MCMCout.indexsave

    if indOut[1]>MCMCout.nsamplesave error("indOut>nsamplesave") end

    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    if Likelihood.SaveMissing == false
        for k in 1:kmax
            MCMCout.mcmc_out["mu"][indOut[1],:,k] = Likelihood.mu.parameteracc[k]
            MCMCout.mcmc_out["nu"][indOut[1],:,k] = Likelihood.nu.parameteracc[k]
            MCMCout.mcmc_out["sigma"][indOut[1],:,k] = vec(Likelihood.sigma.parameteracc[k].mat)
            MCMCout.mcmc_out["pi"][indOut[1],:,k] = Clusterization.pi.parameteracc[k]
            MCMCout.mcmc_out["zeta"][indOut[1],:] = Clusterization.clusterization.zeta[1:(nt-1)]
            MCMCout.mcmc_out["eta"][indOut[1],:,k] = Likelihood.eta.parameteracc[k]
            MCMCout.mcmc_out["rho"][indOut[1],:,k] = Likelihood.rho.parameteracc[k]
        end
    else
        for k in 1:kmax
            MCMCout.mcmc_out["mu"][indOut[1],:,k] = Likelihood.mu.parameteracc[k]
            MCMCout.mcmc_out["nu"][indOut[1],:,k] = Likelihood.nu.parameteracc[k]
            MCMCout.mcmc_out["sigma"][indOut[1],:,k] = vec(Likelihood.sigma.parameteracc[k].mat)
            MCMCout.mcmc_out["pi"][indOut[1],:,k] = Clusterization.pi.parameteracc[k]
            MCMCout.mcmc_out["zeta"][indOut[1],:] = Clusterization.clusterization.zeta[1:(nt-1)]
            MCMCout.mcmc_out["eta"][indOut[1],:,k] = Likelihood.eta.parameteracc[k]
            MCMCout.mcmc_out["rho"][indOut[1],:,k] = Likelihood.rho.parameteracc[k]

            MCMCout.mcmc_out["missing"][indOut[1],:,:]  =  Matrix(transpose(hcat(Likelihood.data.data[Likelihood.miss.indexrow,:]...)))
            # for irow in 1:size(Likelihood.miss.indexrow,1)
            #
            # end
        end
    end

    MCMCout.indexsave[1] += one(MCMCout.indexsave[1])
#Likelihood.data.data[Likelihood.miss.indexrow,:]

#Matrix(Likelihood.data.data[Likelihood.miss.indexrow,:])
    return nothing
end
