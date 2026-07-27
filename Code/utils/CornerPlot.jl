using SpectralFitting
using BSON: @load
using ProgressBars

include("PlottingDefaults.jl")
include("FittingUtils.jl")

const paramLabels = Dict(
    "K" => ["K", "E", "a", "h", "θ", "K2", "a2"], 
    "J" => ["K", "E", "a", "h", "θ", "α13", "ϵ3", "K2", "a2"]
)

function ChiSquaredGrid(result, idx1::Int, idx2::Int, param1, param2)

    stats = Array{Real}(undef, length(param1), length(param2))

    values = copy(result.u)

    pbar = ProgressBar(total = length(param1) * length(param2))

    for (j, p1) in enumerate(param1)
        for (i, p2) in enumerate(param2)

            update(pbar)

            values[idx1] = p1
            values[idx2] = p2

            stats[i,j] = measure(ChiSquared(), result, values)

        end
    end

    stats

end

function SingleCombo(result, idx1::Int, idx2::Int, range1::Real, range2::Real; num=50)

    r1 = result.u[idx1]
    r2 = result.u[idx2]

    param1 = range(r1-range1, r1+range1, num)
    param2 = range(r2-range2, r2+range2, num)

    grid = ChiSquaredGrid(result, idx1, idx2, param1, param2)

    # Finding the standard deviation to plot 1σ,2σ,3σ contours
    stdev = std(grid)
    levels = [stdev, 2*stdev, 3*stdev]

    scatter((r1, r2))
    contour(param1, param2, grid, levels=levels)

end

@load "output/TableJohannsenResultnu80402315002.bson" johannsenResult

plt = SingleCombo(johannsenResult[1], 1, 2, 0.05, 0.05)
display(plt)