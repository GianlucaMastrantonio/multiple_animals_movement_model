
abstract type AbstractMu0 end


struct Mu0DoNotUpdate <: AbstractMu0

end

#### DO UPDATE
abstract type  Mu0DoUpdate  <: AbstractMu0 end

struct Mu0DoUpdate_MvNormal{T<:AbstractMvNormal} <: Mu0DoUpdate

    prior::T
    parameteracc::Vector{Vector{Float64}}
    parameterprop::Vector{Vector{Float64}}
    invmat::AbstractPDMat

    Mu0DoUpdate_MvNormal{T}(prior::T,parameteracc::Vector{Vector{Float64}},parameterprop::Vector{Vector{Float64}},invmat::AbstractPDMat) where {T<:AbstractMvNormal} = new{T}(prior,parameteracc,parameterprop,invmat)
end

function Mu0DoUpdate_MvNormal(params::Tuple)

    prior           = MvNormal(params[1],params[2])
    parameteracc    = params[3]
    parameterprop   = deepcopy(parameteracc)
    invmat          = inv(prior.Σ)

    Mu0DoUpdate_MvNormal{MvNormal}(prior,parameteracc,parameterprop,invmat)

end

#### #### #### #### #### #### #### #### #### ####
#### #### #### sample mu0
#### #### #### #### #### #### #### #### #### ####
# struct Mu0DoUpdate{T<:ContinuousMultivariateDistribution} <: AbstractMu0
#
#     prior::T
#     parameteracc::Vector{Vector{Float64}}
#     parameterprop::Vector{Vector{Float64}}
#
#     Mu0DoUpdate{T}(prior::T,parameteracc::Vector{Vector{Float64}},parameterprop::Vector{Vector{Float64}}) where {T<:ContinuousMultivariateDistribution} = new(prior,parameteracc,parameterprop)
# end
#
# function Mu0DoUpdate_MvNormal(params::Tuple)
#
#     prior           = MvNormal(params[1],params[2])
#     parameteracc    = params[3]
#     parameterprop   = deepcopy(parameteracc)
#     Mu0DoUpdate{typeof(prior)}(prior,parameteracc,parameterprop)
#
# end
