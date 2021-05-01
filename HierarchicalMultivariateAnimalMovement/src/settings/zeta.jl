
abstract type AbstractZeta end


struct ZetaDoNotUpdate <: AbstractZeta
    kmax::Int16

    emptyC::Vector{Int16}
    nonemptyC::Vector{Int16}

    n_obsC::Vector{Int16}

    n_emptyC::Vector{Int16}
    n_nonemptyC::Vector{Int16}

    zeta::Vector{Int16}

    n_itojC::Vector{Vector{Int16}}

    ZetaDoNotUpdate(kmax::Int16,emptyC::Vector{Int16},nonemptyC::Vector{Int16},n_obsC::Vector{Int16},n_emptyC::Vector{Int16},n_nonemptyC::Vector{Int16},zeta::Vector{Int16},n_itojC::Vector{Vector{Int16}}) = new(kmax,emptyC,nonemptyC,n_obsC,n_emptyC,n_nonemptyC,zeta,n_itojC)
end

function ZetaDoNotUpdate(zeta::Vector{Int16},kmax::Int16)

    if maximum(zeta)>kmax error("max of zeta must be less or equal than kmax") end
    if minimum(zeta)<1 error("min of zeta must be greater or equal to than 1") end
    #
    # emptyC      = Vector{Int16}(undef,kmax)
    # nonemptyC   = Vector{Int16}(undef,kmax)
    #
    # n_obsC       = Vector{Int16}(undef,kmax)
    #
    # n_itojC     = Vector{Vector{Int16}}()
    #

    emptyC      = zeros(Int16,kmax)
    nonemptyC   = zeros(Int16,kmax)

    n_obsC       = zeros(Int16,kmax)

    n_itojC     = Vector{Vector{Int16}}()

    for k in 1:kmax
        push!(n_itojC, zeros(Int16,kmax))
    end

    n_nonemptyC_1,n_emptyC_1 = Update_ObsInClust!(kmax,zeta,emptyC,nonemptyC,n_obsC,n_itojC)

    n_nonemptyC = Int16[n_nonemptyC_1]
    n_emptyC    =  Int16[n_emptyC_1]
    ZetaDoNotUpdate(kmax,emptyC,nonemptyC,n_obsC,n_emptyC,n_nonemptyC,zeta,n_itojC)

end

#### DO UPDATE
struct ZetaDoUpdate <: AbstractZeta

    kmax::Int16

    emptyC::Vector{Int16}
    nonemptyC::Vector{Int16}

    n_obsC::Vector{Int16}

    n_emptyC::Vector{Int16}
    n_nonemptyC::Vector{Int16}

    zeta::Vector{Int16}

    n_itojC::Vector{Vector{Int16}}

    ZetaDoUpdate(kmax::Int16,emptyC::Vector{Int16},nonemptyC::Vector{Int16},n_obsC::Vector{Int16},n_emptyC::Vector{Int16},n_nonemptyC::Vector{Int16},zeta::Vector{Int16},n_itojC::Vector{Vector{Int16}}) = new(kmax,emptyC,nonemptyC,n_obsC,n_emptyC,n_nonemptyC,zeta,n_itojC)
end

function ZetaDoUpdate(zeta::Vector{Int16},kmax::Int16)

    if maximum(zeta)>kmax error("max of zeta must be less or equal than kmax") end
    if minimum(zeta)<1 error("min of zeta must be greater or equal to than 1") end

    emptyC      = zeros(Int16,kmax)
    nonemptyC   = zeros(Int16,kmax)

    n_obsC       = zeros(Int16,kmax)

    n_itojC     = Vector{Vector{Int16}}()

    for k in 1:kmax
        push!(n_itojC, zeros(Int16,kmax))
    end

    n_nonemptyC_1,n_emptyC_1 = Update_ObsInClust!(kmax,zeta,emptyC,nonemptyC,n_obsC,n_itojC)

    n_nonemptyC = Int16[n_nonemptyC_1]
    n_emptyC    =  Int16[n_emptyC_1]
    ZetaDoUpdate(kmax,emptyC,nonemptyC,n_obsC,n_emptyC,n_nonemptyC,zeta,n_itojC)

end

function Update_ObsInClust!(kmax::Integer,zeta::Vector{Int16},emptyC::Vector{Int16},nonemptyC::Vector{Int16},n_obsC::Vector{Int16},n_itojC::Vector{Vector{Int16}})::Tuple{Int16,Int16}

    for k in 1:kmax
        @inbounds emptyC[k]       = Int16(0)
        @inbounds nonemptyC[k]    = Int16(0)
        @inbounds n_obsC[k]       = Int16(0)

        for k1 in 1:kmax
            @inbounds n_itojC[k][k1] = Int16(0)
        end

    end

    zprev = 1
    for i in 1:(size(zeta,1)-1)
        @inbounds n_itojC[zprev][zeta[i]] += one(zeta[i])
        @inbounds n_obsC[zeta[i]]         += one(zeta[i])
        @inbounds zprev                   = zeta[i]
    end
    # i = nt
    # n_itojC[zprev][zeta[i]] += one(zeta[i])

    iempty          = Int16(0)
    inonempty       = Int16(0)
    for k in 1:kmax
        if n_obsC[k] != 0
            @inbounds inonempty += one(inonempty)
            @inbounds nonemptyC[inonempty] = k
        else
            @inbounds iempty += one(iempty)
            @inbounds emptyC[iempty] = k
        end
    end
    return (inonempty , iempty)
end

function Update_ObsInClust!(zetastruct::ZetaDoUpdate)
    n_nonemptyC,n_emptyC = Update_ObsInClust!(zetastruct.kmax,zetastruct.zeta,zetastruct.emptyC,zetastruct.nonemptyC,zetastruct.n_obsC,zetastruct.n_itojC)

    zetastruct.n_emptyC[1,1]     = n_emptyC
    zetastruct.n_nonemptyC[1,1]  = n_nonemptyC
end
