
abstract type AbstractPosDefMatPar end
abstract type AbstractPosDefMatParDoUpdate <:AbstractPosDefMatPar  end
struct        PosDefMatParDoNotUpdate      <:AbstractPosDefMatPar   end

#### #### #### #### #### #### #### #### #### ####
#### DO UPDATE
#### #### #### #### #### #### #### #### #### ####
struct PosDefMatInverseWishart <: AbstractPosDefMatParDoUpdate

    prior::InverseWishart{Float64,PDMat{Float64,Array{Float64,2}}}
    parameteracc::Vector{PDMat{Float64,Array{Float64,2}}}
    parameterprop::Vector{PDMat{Float64,Array{Float64,2}}}

    parameteraccinv::Vector{PDMat{Float64,Array{Float64,2}}}
    parameterpropinv::Vector{PDMat{Float64,Array{Float64,2}}}

    PosDefMatInverseWishart(
    prior::InverseWishart{Float64,PDMat{Float64,Array{Float64,2}}},
    parameteracc::Vector{PDMat{Float64,Array{Float64,2}}},
    parameterprop::Vector{PDMat{Float64,Array{Float64,2}}},

    parameteraccinv::Vector{PDMat{Float64,Array{Float64,2}}},
    parameterpropinv::Vector{PDMat{Float64,Array{Float64,2}}}

    )  = new(prior,parameteracc,parameterprop,parameteraccinv,parameterpropinv)
end

function PosDefMatInverseWishart(params::Tuple{Float64,Array{Float64,2},Array{Array{Float64,2},1}})

    prior           = InverseWishart(params[1],params[2])

    parameteracc    = Vector{PDMat{Float64,Array{Float64,2}}}()
    for k in 1:size(params[3],1)
        push!(parameteracc,PDMat(params[3][k]))
    end
    parameterprop   = deepcopy(parameteracc)

    parameteraccinv = Vector{PDMat{Float64,Array{Float64,2}}}()
    for k in 1:size(params[3],1)
        push!(parameteraccinv,inv(parameteracc[k]))
    end
    parameterpropinv   = deepcopy(parameteraccinv)

    PosDefMatInverseWishart(prior,parameteracc,parameterprop,parameteraccinv,parameterpropinv)

end
