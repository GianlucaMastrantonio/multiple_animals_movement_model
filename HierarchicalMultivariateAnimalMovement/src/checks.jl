

### mu prior parameters and inits for mixture
function check_MixtureMvNormal(name::String,inits::Vector{Vector{Float64}}, prior::Dict, kmax::Integer, dimpar::Integer)

    # parameters
    μ = prior["μ"]::Vector{Float64}
    Σ = prior["Σ"]::Matrix{Float64}

    # mu
    if !haskey(prior, "μ") error(string("μ not in ", name, "prior"))  end
    dimμ_1    = size(μ,1)
    dimμ_2    = size(μ,2)
    Distributions.@check_args(name, dimμ_1 == dimpar)

    # Sigma
    if !haskey(prior, "Σ") error(string("Σ not in ", name, "prior"))  end
    dimΣ_1    = size(Σ,1)
    dimΣ_2    = size(Σ,2)
    Distributions.@check_args(name, dimΣ_1 == dimpar)
    Distributions.@check_args(name, dimΣ_2 == dimpar)
    Distributions.@check_args(name, isposdef(Σ))

    # inits
    diminits_1  = size(inits,1)
    Distributions.@check_args(name, diminits_1 == kmax)
    for i in 1:diminits_1
        diminits_2  = size(inits[i],1)
        Distributions.@check_args(name, diminits_2 == dimpar)
    end


    return (μ,Σ,inits)
end


function check_MixtureMvBeta(name::String,inits::Vector{Vector{Float64}}, prior::Dict, kmax::Integer, dimpar::Integer)

    # parameters
    α = prior["α"]::Vector{Float64}
    β = prior["β"]::Vector{Float64}

    # alpha
    if !haskey(prior, "α") error(string("α not in ", name, "prior"))  end
    dimα_1    = size(α,1)
    dimα_2    = size(α,2)
    Distributions.@check_args(name, dimα_1 == dimpar)

    # beta
    if !haskey(prior, "β") error(string("β not in ", name, "prior"))  end
    dimβ_1    = size(β,1)
    dimβ_2    = size(β,2)
    Distributions.@check_args(name, dimβ_1 == dimpar)

    # inits
    diminits_1  = size(inits,1)
    Distributions.@check_args(name, diminits_1 == kmax)
    for i in 1:diminits_1
        diminits_2  = size(inits[i],1)
        Distributions.@check_args(name, diminits_2 == dimpar)
    end


    return (α,β,inits)
end




function check_MixtureMvUniform(name::String,inits::Vector{Vector{Float64}}, prior::Dict, kmax::Integer, dimpar::Integer)

    # parameters
    a = prior["a"]::Vector{Float64}
    b = prior["b"]::Vector{Float64}

    # alpha
    if !haskey(prior, "a") error(string("a not in ", name, "prior"))  end
    dima_1    = size(a,1)
    dima_2    = size(a,2)
    Distributions.@check_args(name, dima_1 == dimpar)

    # beta
    if !haskey(prior, "b") error(string("b not in ", name, "prior"))  end
    dimb_1    = size(b,1)
    dimb_2    = size(b,2)
    Distributions.@check_args(name, dimb_1 == dimpar)

    # inits
    diminits_1  = size(inits,1)
    Distributions.@check_args(name, diminits_1 == kmax)
    for i in 1:diminits_1
        diminits_2  = size(inits[i],1)
        Distributions.@check_args(name, diminits_2 == dimpar)
    end


    return (a,b,inits)
end



function check_MixtureInverseWishart(name::String,inits::Vector{Matrix{Float64}}, prior::Dict, kmax::Integer, dimpar::Integer)

    # parameters

    df = prior["df"]::Float64
    Ψ  = prior["Ψ"]::Matrix{Float64}

    # Ψ
    if !haskey(prior, "Ψ") error(string("Ψ not in ", name, "prior"))  end
    dimΨ_1    = size(Ψ,1)
    dimΨ_2    = size(Ψ,2)
    Distributions.@check_args(name, dimΨ_1 == dimpar)
    Distributions.@check_args(name, dimΨ_2 == dimpar)


    # inits
    diminits_1  = size(inits,1)
    Distributions.@check_args(name, diminits_1 == kmax)
    for i in 1:diminits_1
        diminits_2  = size(inits[i],1)
        diminits_3  = size(inits[i],2)
        Distributions.@check_args(name, diminits_2 == dimpar)
        Distributions.@check_args(name, diminits_3 == dimpar)
    end


    return (df,Ψ,inits)
end

function check_ClusterizationDirichlet(name::String,inits::Vector{Vector{Float64}}, prior::Dict, kmax::Integer)

    # parameters

    alpha = prior["alpha"]::Vector{Float64}


    # alpha
    if !haskey(prior, "alpha") error(string("alpha not in ", name, "prior"))  end
    dimalpha_1    = size(alpha,1)
    Distributions.@check_args(name, dimalpha_1 == kmax)

    # inits
    diminits_1  = size(inits,1)
    Distributions.@check_args(name, diminits_1 == kmax)
    for i in 1:diminits_1
        diminits_2  = size(inits[i],1)
        Distributions.@check_args(name, diminits_2 == kmax)
    end

    return (alpha,inits)
end
