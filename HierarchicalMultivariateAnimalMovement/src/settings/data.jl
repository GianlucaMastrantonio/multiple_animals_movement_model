

#### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### DATASET
#### #### #### #### #### #### #### #### #### #### #### ####
abstract type AbstractDataset end

struct CoordinatesDataset <: AbstractDataset

    data::Vector{Vector{Union{Float64}}}
    nt::Int32
    ncol::Int32
    nanimals::Int16

    CoordinatesDataset(data::Vector{Vector{Union{Float64}}},nt::Int32,ncol::Int32,nanimals::Int16) = new(data,nt,ncol,nanimals)
end

function CoordinatesDataset(datacopy::Matrix{Union{Missing,Float64}})

    nt          = Int32(size(datacopy)[1])
    ncol        = Int32(size(datacopy)[2])
    nanimals    = Int16((ncol/2))
    data        = Vector{Vector{Union{Float64}}}()
    for i in 1:nt
        push!(data, zeros(Union{Float64},ncol))
        for j in 1:ncol
            data[i][j] = datacopy[i,j]
        end
    end


    CoordinatesDataset(data,nt,ncol,nanimals)
end

#### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### MODEL
#### #### #### #### #### #### #### #### #### #### ####


####
abstract type AbstractLikelihood end

struct Likelihood_OUmodel{
            Tmiss<:AbstractMissing,
            Tmu<:AbstractVecPar,
            Tnu<:AbstractVecPar,
            Tsigma<:AbstractPosDefMatPar ,
            Tzeta<:AbstractZeta
                        }     <: AbstractLikelihood

    data::CoordinatesDataset
    miss::Tmiss
    mu::Tmu
    nu::Tnu
    sigma::Tsigma
    clusterization::Tzeta
    kmax::Int16

    Likelihood_OUmodel{Tmiss,Tmu,Tnu,Tsigma,Tzeta}(data::CoordinatesDataset,miss::Tmiss,mu::Tmu,nu::Tnu,sigma::Tsigma,clusterization::Tzeta,kmax::Int16) where {
                Tmiss<:AbstractMissing,
                Tmu<:AbstractVecPar,
                Tnu<:AbstractVecPar,
                Tsigma<:AbstractPosDefMatPar,
                Tzeta<:AbstractZeta
                            } = new{Tmiss,Tmu,Tnu,Tsigma,Tzeta}(data,miss,mu,nu,sigma,clusterization,kmax)
end

function Likelihood_OUmodel(data::CoordinatesDataset, miss::Tmiss,mu::Tmu,nu::Tnu,sigma::Tsigma,clusterization::Tzeta,kmax::Int16) where {
            Tmiss<:AbstractMissing,
            Tmu<:AbstractVecPar,
            Tnu<:AbstractVecPar,
            Tsigma<:AbstractPosDefMatPar,
            Tzeta<:AbstractZeta
                        }
    Likelihood_OUmodel{Tmiss,Tmu,Tnu,Tsigma,Tzeta}(data,miss,mu,nu,sigma,clusterization,kmax)

end

####
struct Likelihood_OU_CircLinmodel{
            Tmiss<:AbstractMissing,
            Tmu<:AbstractVecPar,
            Tnu<:AbstractVecPar,
            Tsigma<:AbstractPosDefMatPar ,
            Tzeta<:AbstractZeta,
            Teta<:AbstractVecPar,
            Trho<:AbstractVecPar,
            Tangle<:AbstractVecPar
                        }     <: AbstractLikelihood

    data::CoordinatesDataset
    miss::Tmiss
    mu::Tmu
    nu::Tnu
    sigma::Tsigma
    clusterization::Tzeta
    kmax::Int16
    eta::Teta
    rho::Trho
    Angle::Tangle
    SaveMissing::Bool

    Likelihood_OU_CircLinmodel{Tmiss,Tmu,Tnu,Tsigma,Tzeta,Teta,Trho,Tangle}(data::CoordinatesDataset,miss::Tmiss,mu::Tmu,nu::Tnu,sigma::Tsigma,clusterization::Tzeta,kmax::Int16,eta::Teta,rho::Trho,Angle::Tangle,SaveMissing::Bool) where {
                Tmiss<:AbstractMissing,
                Tmu<:AbstractVecPar,
                Tnu<:AbstractVecPar,
                Tsigma<:AbstractPosDefMatPar,
                Tzeta<:AbstractZeta,
                Teta<:AbstractVecPar,
                Trho<:AbstractVecPar,
                Tangle<:AbstractVecPar
                            } = new{Tmiss,Tmu,Tnu,Tsigma,Tzeta,Teta,Trho,Tangle}(data,miss,mu,nu,sigma,clusterization,kmax,eta,rho,Angle,SaveMissing)
end



function Likelihood_OU_CircLinmodel(data::CoordinatesDataset, miss::Tmiss,mu::Tmu,nu::Tnu,sigma::Tsigma,clusterization::Tzeta,kmax::Int16,eta::Teta,rho::Trho,Angle::Tangle,SaveMissing::Bool) where {
            Tmiss<:AbstractMissing,
            Tmu<:AbstractVecPar,
            Tnu<:AbstractVecPar,
            Tsigma<:AbstractPosDefMatPar,
            Tzeta<:AbstractZeta,
            Teta<:AbstractVecPar,
            Trho<:AbstractVecPar,
            Tangle<:AbstractVecPar
                        }
    Likelihood_OU_CircLinmodel{Tmiss,Tmu,Tnu,Tsigma,Tzeta,Teta,Trho,Tangle}(data,miss,mu,nu,sigma,clusterization,kmax,eta,rho,Angle,SaveMissing)

end


#### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #SECT HMM models
#### #### #### #### #### #### #### #### #### #### ####
#
# abstract type HMMtype end
# struct Seperated <: HMMtype end
# struct Joint     <: HMMtype end
#
# abstract type DPtype end
# struct Simple        <: DPtype end
# struct Hierarichical <: DPtype end
#
# abstract type ModelSecondLevel{F<:HMMtype,S<:DPtype} end
#
# const HMM_sep_simp     = ModelSecondLevel{Seperated,Simple}
# const HMM_sep_hierc    = ModelSecondLevel{Seperated,Hierarichical}
# const HMM_joint_simp   = ModelSecondLevel{Joint,Simple}
# const HMM_joint_hierc  = ModelSecondLevel{Joint,Hierarichical}
#
# const HMM_sep{S<:DPtype}          = ModelSecondLevel{Seperated,S}
# const HMM_joint{S<:DPtype}        = ModelSecondLevel{Joint,S}
# const HMM_simple{F<:HMMtype}      = ModelSecondLevel{F,Simple}
# const HMM_hierc{F<:HMMtype}       = ModelSecondLevel{F,Hierarichical}
abstract type AbstractClusterization end

struct Clusterization_HMM{Tpi<:AbstractVecPar,Tzeta<:AbstractZeta}   <: AbstractClusterization

    clusterization::Tzeta
    pi::Tpi
    mcmc_initpi::Vector{Float64}
    mcmc_initz::Vector{Int16}

    Clusterization_HMM{Tpi,Tzeta}(clusterization::Tzeta,pi::Tpi,mcmc_initpi,mcmc_initz) where {Tpi<:AbstractVecPar,Tzeta<:AbstractZeta} = new{Tpi,Tzeta}(clusterization,pi,mcmc_initpi,mcmc_initz)
end

function Clusterization_HMM(clusterization::Tzeta,pi::Tpi,mcmc_initpi::Vector{Float64},mcmc_initz::Vector{Int16}) where {Tpi<:AbstractVecPar,Tzeta<:AbstractZeta}
    Clusterization_HMM{Tpi,Tzeta}(clusterization,pi,mcmc_initpi,mcmc_initz)
end





###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ######

abstract type AbstractClusterization_Divided  <: AbstractClusterization end
abstract type AbstractClusterization_Joint  <: AbstractClusterization end


struct Clusterization_HDPHMM_Joint{Tpi<:AbstractVecPar,Tzeta<:AbstractZeta}   <: AbstractClusterization_Joint

    clusterization::Tzeta
    pi::Tpi
    mcmc_initpi::Vector{Float64}
    mcmc_initz::Vector{Int16}
    mcmc_beta::Vector{Float64}
    mcmc_rho::Vector{Float64}
    mcmc_gamma::Vector{Float64}
    mcmc_ak::Vector{Float64}

    prior_ak::Vector{Float64}
    prior_gamma::Vector{Float64}
    prior_rho::Vector{Float64}
    MatM::Matrix{Float64}
    MatMbar::Matrix{Float64}

    Clusterization_HDPHMM_Joint{Tpi,Tzeta}(clusterization::Tzeta,pi::Tpi,mcmc_initpi,mcmc_initz,mcmc_beta::Vector{Float64},mcmc_rho::Vector{Float64},mcmc_gamma::Vector{Float64},  mcmc_ak::Vector{Float64},prior_ak::Vector{Float64},prior_gamma::Vector{Float64},prior_rho::Vector{Float64},MatM::Matrix{Float64},MatMbar::Matrix{Float64}) where {Tpi<:AbstractVecPar,Tzeta<:AbstractZeta} = new{Tpi,Tzeta}(clusterization,pi,mcmc_initpi,mcmc_initz,mcmc_beta,mcmc_rho,mcmc_gamma,  mcmc_ak,prior_ak,prior_gamma,prior_rho,MatM,MatMbar)
end

function Clusterization_HDPHMM_Joint(clusterization::Tzeta, pi::Tpi,mcmc_initpi::Vector{Float64},mcmc_initz::Vector{Int16},mcmc_beta::Vector{Float64},mcmc_rho::Vector{Float64},mcmc_gamma::Vector{Float64},  mcmc_ak::Vector{Float64},prior_ak::Vector{Float64},prior_gamma::Vector{Float64},prior_rho::Vector{Float64}) where {Tpi<:AbstractVecPar,Tzeta<:AbstractZeta}

    ktot        = size(mcmc_beta)[1]
    #print(ktot)
    MatM        = zeros(Float64,ktot,ktot)
    MatMbar     = zeros(Float64,ktot,ktot)

    Clusterization_HDPHMM_Joint{Tpi,Tzeta}(clusterization,pi,mcmc_initpi,mcmc_initz,mcmc_beta,mcmc_rho,mcmc_gamma,  mcmc_ak,prior_ak,prior_gamma,prior_rho,MatM,MatMbar)
end


struct Clusterization_HDPHMM{Tpi<:AbstractVecPar,Tzeta<:AbstractZeta}   <: AbstractClusterization_Divided

    clusterization::Tzeta
    pi::Tpi
    mcmc_initpi::Vector{Float64}
    mcmc_initz::Vector{Int16}
    mcmc_beta::Vector{Float64}
    mcmc_rho::Vector{Float64}
    mcmc_gamma::Vector{Float64}
    mcmc_ak::Vector{Float64}

    prior_ak::Vector{Float64}
    prior_gamma::Vector{Float64}
    prior_rho::Vector{Float64}
    MatM::Matrix{Float64}
    MatMbar::Matrix{Float64}
    gamma_par::Vector{Float64}

    Clusterization_HDPHMM{Tpi,Tzeta}(clusterization::Tzeta,pi::Tpi,mcmc_initpi,mcmc_initz,mcmc_beta::Vector{Float64},mcmc_rho::Vector{Float64},mcmc_gamma::Vector{Float64},  mcmc_ak::Vector{Float64},prior_ak::Vector{Float64},prior_gamma::Vector{Float64},prior_rho::Vector{Float64},MatM::Matrix{Float64},MatMbar::Matrix{Float64},gamma_par::Vector{Float64}) where {Tpi<:AbstractVecPar,Tzeta<:AbstractZeta} = new{Tpi,Tzeta}(clusterization,pi,mcmc_initpi,mcmc_initz,mcmc_beta,mcmc_rho,mcmc_gamma,  mcmc_ak,prior_ak,prior_gamma,prior_rho,MatM,MatMbar,gamma_par)
end

function Clusterization_HDPHMM(clusterization::Tzeta, pi::Tpi,mcmc_initpi::Vector{Float64},mcmc_initz::Vector{Int16},mcmc_beta::Vector{Float64},mcmc_rho::Vector{Float64},mcmc_gamma::Vector{Float64},  mcmc_ak::Vector{Float64},prior_ak::Vector{Float64},prior_gamma::Vector{Float64},prior_rho::Vector{Float64}) where {Tpi<:AbstractVecPar,Tzeta<:AbstractZeta}

    ktot        = size(mcmc_beta)[1]
    #print(ktot)
    MatM        = zeros(Float64,ktot,ktot)
    MatMbar     = zeros(Float64,ktot,ktot)
    gamma_par   = ones(Float64,5)

    Clusterization_HDPHMM{Tpi,Tzeta}(clusterization,pi,mcmc_initpi,mcmc_initz,mcmc_beta,mcmc_rho,mcmc_gamma,  mcmc_ak,prior_ak,prior_gamma,prior_rho,MatM,MatMbar,gamma_par)
end



#### #### #### #### #### #### #### #### #### #### ####
#### ####
#### #### #### #### #### #### #### #### #### #### ####

abstract type Parameter_Hierarchical end

struct PosDefMatParParameter_Hierarchical <:Parameter_Hierarchical

    clust::Matrix{Int16}
    prob::Vector{Float64}
    par::AbstractPosDefMatPar
    obsinclust::Matrix{Int16}

    PosDefMatParParameter_Hierarchical(clust::Matrix{Int16}, prob::Vector{Float64},par::AbstractPosDefMatPar, obsinclust::Matrix{Int16}) = new(clust, prob,par, obsinclust)
end

function PosDefMatParParameter_Hierarchical(par_app::AbstractPosDefMatPar,nanim::Int16)

    kmax        = size(par_app.parameteracc)[1]
    clust       = Matrix{Int16}(undef,kmax,nanim )
    for i in 1:nanim
        clust[:,i] = rem.((1:kmax) .-1   .+ (i-1)*20 ,kmax) .+1
    end
    obsinclust  = zeros(Int16,kmax,nanim)
    prob		= ones(Float64,kmax)./kmax
    par         = deepcopy(par_app)

	PosDefMatParParameter_Hierarchical(clust, prob,par, obsinclust)
end

struct VecParParameter_Hierarchical <:Parameter_Hierarchical

    clust::Matrix{Int16}
    prob::Vector{Float64}
    par::AbstractVecPar
    obsinclust::Matrix{Int16}

    VecParParameter_Hierarchical(clust::Matrix{Int16}, prob::Vector{Float64},par::AbstractVecPar, obsinclust::Matrix{Int16}) = new(clust, prob,par, obsinclust)
end

function VecParParameter_Hierarchical(par_app::AbstractVecPar, nanim::Int16)

	kmax        = size(par_app.parameteracc)[1]
    clust       = Matrix{Int16}(undef,kmax,nanim )
    for i in 1:nanim
        clust[:,i] = rem.((1:kmax) .-1   .+ (i-1)*20 ,kmax) .+1
    end
    obsinclust  = zeros(Int16,kmax,nanim)
    prob		= ones(Float64,kmax)./kmax
    par         = deepcopy(par_app)

	VecParParameter_Hierarchical(clust, prob,par, obsinclust)
end


struct OptionHierarchicalParameters

    h_mu::VecParParameter_Hierarchical
	h_eta::VecParParameter_Hierarchical
	h_nu::VecParParameter_Hierarchical
	h_rho::VecParParameter_Hierarchical
	h_sigma::PosDefMatParParameter_Hierarchical

    OptionHierarchicalParameters(h_mu::VecParParameter_Hierarchical, h_eta::VecParParameter_Hierarchical, h_nu::VecParParameter_Hierarchical, h_rho::VecParParameter_Hierarchical, h_sigma::PosDefMatParParameter_Hierarchical) = new(h_mu, h_eta, h_nu, h_rho, h_sigma)

end


function OptionHierarchicalParameters(likelihood::Vector{Likelihood_OU_CircLinmodel})

	nanim = size(likelihood)[1]

	h_mu    = VecParParameter_Hierarchical(likelihood[1].mu, Int16(nanim))
	h_eta   = VecParParameter_Hierarchical(likelihood[1].eta, Int16(nanim))
	h_nu    = VecParParameter_Hierarchical(likelihood[1].nu, Int16(nanim))
	h_rho   = VecParParameter_Hierarchical(likelihood[1].rho, Int16(nanim))
	h_sigma = PosDefMatParParameter_Hierarchical(likelihood[1].sigma, Int16(nanim))

	OptionHierarchicalParameters(h_mu, h_eta, h_nu, h_rho, h_sigma)

end
