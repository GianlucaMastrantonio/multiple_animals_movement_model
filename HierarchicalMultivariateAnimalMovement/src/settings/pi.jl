
abstract type AbstractPi end
struct PiDoNotUpdate <: AbstractPi

end

#### DO UPDATE
struct PiDoUpdate{T<:ContinuousMultivariateDistribution} <: AbstractPi

    prior::T
    parameteracc::Vector{Vector{Float64}}
    parameterprop::Vector{Vector{Float64}}

    PiDoUpdate{T}(prior::T,parameteracc::Vector{Vector{Float64}},parameterprop::Vector{Vector{Float64}}) where {T<:ContinuousMultivariateDistribution} = new(prior,parameteracc,parameterprop)
end

function PiDoUpdate_Dirichlet(params::Tuple)

    prior           = Dirichlet(params[1])
    parameteracc    = params[2]
    parameterprop   = deepcopy(parameteracc)
    PiDoUpdate{typeof(prior)}(prior,parameteracc,parameterprop)

end
