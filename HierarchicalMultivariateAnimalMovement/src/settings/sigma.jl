
abstract type AbstractSigma end


struct SigmaDoNotUpdate <: AbstractSigma

end

#### DO UPDATE
struct SigmaDoUpdate{T<:ContinuousMatrixDistribution} <: AbstractSigma

    prior::T
    parameteracc::Vector{Matrix{Float64}}
    parameterprop::Vector{Matrix{Float64}}

    SigmaDoUpdate{T}(prior::T,parameteracc::Vector{Matrix{Float64}},parameterprop::Vector{Matrix{Float64}}) where {T<:ContinuousMatrixDistribution} = new(prior,parameteracc,parameterprop)
end

function SigmaDoUpdate_InverseWishart(params::Tuple)

    prior   = InverseWishart(params[1],params[2])
    parameteracc    = params[3]
    parameterprop   = deepcopy(parameteracc)

    SigmaDoUpdate{typeof(prior)}(prior,parameteracc,parameterprop)

end
