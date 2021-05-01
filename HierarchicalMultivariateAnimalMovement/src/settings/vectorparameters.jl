
abstract type   AbstractVecPar                                     end
abstract type   AbstractVecParDoUpdate <:AbstractVecPar   end
struct          VecParDoNotUpdate      <:AbstractVecPar

    #parameteracc::Vector{Vector{Float64}}
    #parameterprop::Vector{Vector{Float64}}

end

# function VecParDoNotUpdate(params::Vector{Vector{Float64}})
#
#
#     parameteracc    = params
#     parameterprop   = deepcopy(parameteracc)
#
#     VecParDoNotUpdate(parameteracc,parameterprop)
#
# end

#### #### #### #### #### #### #### #### #### ####
#### DO UPDATE
#### #### #### #### #### #### #### #### #### ####
##### MULTI UNIFORM
struct VecParMvNoPrior <: AbstractVecParDoUpdate

    #prior::Product{Continuous,Uniform{Float64},Array{Uniform{Float64},1}}
    parameteracc::Vector{Vector{Float64}}
    parameterprop::Vector{Vector{Float64}}

    VecParMvNoPrior(parameteracc::Vector{Vector{Float64}},parameterprop::Vector{Vector{Float64}}) = new(parameteracc,parameterprop)
end

function VecParMvNoPrior(params::Vector{Vector{Float64}})


    parameteracc    = params
    parameterprop   = deepcopy(parameteracc)

    VecParMvNoPrior(parameteracc,parameterprop)

end


##### MULTI NORMAL
abstract type   AbstractVecParMvNormal <:  AbstractVecParDoUpdate end
struct VecParMvNormal <: AbstractVecParMvNormal

    prior::MvNormal{Float64,PDMat{Float64,Array{Float64,2}},Array{Float64,1}}
    parameteracc::Vector{Vector{Float64}}
    parameterprop::Vector{Vector{Float64}}
    prior_invmat::PDMat{Float64,Array{Float64,2}}

    VecParMvNormal(prior::MvNormal{Float64,PDMat{Float64,Array{Float64,2}},Array{Float64,1}}, parameteracc::Vector{Vector{Float64}},parameterprop::Vector{Vector{Float64}},prior_invmat::PDMat{Float64,Array{Float64,2}})  = new(prior,parameteracc,parameterprop,prior_invmat)

#VecParMvNormal(::MvNormal{Float64,PDMat{Float64,Array{Float64,2}},Array{Float64,1}}, ::Array{Array{Float64,1},1}, ::Array{Array{Float64,1},1}, ::PDMat{Float64,Array{Float64,2}})
end

function VecParMvNormal(params::Tuple{Vector{Float64},Matrix{Float64},Vector{Vector{Float64}}})

    prior           = MvNormal(params[1],params[2])
    parameteracc    = params[3]
    parameterprop   = deepcopy(parameteracc)
    prior_invmat    = inv(prior.Σ)

    VecParMvNormal(prior,parameteracc,parameterprop,prior_invmat)

end



##### MULTI BETA
struct VecParMvBeta <: AbstractVecParDoUpdate

    prior::Product{Continuous,Beta{Float64},Array{Beta{Float64},1}}
    parameteracc::Vector{Vector{Float64}}
    parameterprop::Vector{Vector{Float64}}

    VecParMvBeta(prior::Product{Continuous,Beta{Float64},Array{Beta{Float64},1}},parameteracc::Vector{Vector{Float64}},parameterprop::Vector{Vector{Float64}}) = new(prior,parameteracc,parameterprop)
end

function VecParMvBeta(params::Tuple{Vector{Float64},Vector{Float64},Vector{Vector{Float64}}})

    ds      = Beta.(params[1],params[2])
    prior   = product_distribution(ds)
    parameteracc    = params[3]
    parameterprop   = deepcopy(parameteracc)

    VecParMvBeta(prior,parameteracc,parameterprop)

end


##### MULTI UNIFORM
struct VecParMvUniform <: AbstractVecParDoUpdate

    prior::Product{Continuous,Uniform{Float64},Array{Uniform{Float64},1}}
    parameteracc::Vector{Vector{Float64}}
    parameterprop::Vector{Vector{Float64}}

    VecParMvUniform(prior::Product{Continuous,Uniform{Float64},Array{Uniform{Float64},1}},parameteracc::Vector{Vector{Float64}},parameterprop::Vector{Vector{Float64}}) = new(prior,parameteracc,parameterprop)
end

function VecParMvUniform(params::Tuple{Vector{Float64},Vector{Float64},Vector{Vector{Float64}}})

    ds      = Uniform.(params[1],params[2])
    prior   = product_distribution(ds)
    parameteracc    = params[3]
    parameterprop   = deepcopy(parameteracc)

    VecParMvUniform(prior,parameteracc,parameterprop)

end


##### DIRICHLET
struct VecParDirichlet <: AbstractVecParDoUpdate

    prior::Dirichlet{Float64}
    parameteracc::Vector{Vector{Float64}}
    parameterprop::Vector{Vector{Float64}}

    VecParDirichlet(prior::Dirichlet{Float64},parameteracc::Vector{Vector{Float64}},parameterprop::Vector{Vector{Float64}}) = new(prior,parameteracc,parameterprop)
end

function VecParDirichlet(params::Tuple{Vector{Float64},Vector{Vector{Float64}}})

    prior           = Dirichlet(params[1])
    parameteracc    = params[2]
    parameterprop   = deepcopy(parameteracc)
    VecParDirichlet(prior,parameteracc,parameterprop)

end
