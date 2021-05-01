abstract type AbstractMissing end


#### DO NOT UPDATE
struct MissingDoNotUpdate<: AbstractMissing
    doupdate::Bool

    function MissingDoNotUpdate()
        doupdate = false
        new(doupdate)
    end
end

function MissingDoNotUpdate(dataset::Matrix{Union{Float64,Missing}})

    indexcol = Vector{Vector{Int32}}(undef,0)
    indexrow = Vector{Int32}(undef,0)

    for irow in 1:size(dataset)[1]
        missapp = findall(ismissing, dataset[irow,:])
        if(size(missapp)[1]!==0)
            push!(indexrow ,irow)
            push!(indexcol,missapp)
        end
    end
    datasetApp = Impute.interp(dataset) |> Impute.locf() |> Impute.nocb()

    for icol in 1:size(dataset,2)
        for irow in 1:size(dataset,1)
            dataset[irow,icol] = datasetApp[irow,icol]
        end
    end
    MissingDoNotUpdate()
end

#### DO UPDATE
struct MissingDoUpdate <: AbstractMissing

    doupdate::Bool
    indexrow::Vector{Int32}
    indexcol::Vector{Vector{Int32}}


    MissingDoUpdate(doupdate::Bool,MissingIndexRow::Vector{Int32},MissingIndexCol::Vector{Vector{Int32}}) = new(doupdate,MissingIndexRow,MissingIndexCol)
end

function MissingDoUpdate(dataset::Matrix{Union{Float64,Missing}})

    indexcol = Vector{Vector{Int32}}(undef,0)
    indexrow = Vector{Int32}(undef,0)

    for irow in 1:size(dataset)[1]
        missapp = findall(ismissing, dataset[irow,:])
        if(size(missapp)[1]!==0)
            push!(indexrow ,irow)
            push!(indexcol,missapp)
        end
    end
    datasetApp = Impute.interp(dataset) |> Impute.locf() |> Impute.nocb()

    for icol in 1:size(dataset,2)
        for irow in 1:size(dataset,1)
            dataset[irow,icol] = datasetApp[irow,icol]
        end
    end
    MissingDoUpdate(true,indexrow ,indexcol)
end
