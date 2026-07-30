using SpectralFitting
using BSON: @load
using ProgressBars
using StatsBase: std
using Plots

pyplot()

cm2in(x) = 0.3937008x
cm2px(x) = Int(round(cm2in(x) * 100))

default(
    # Fonts
    titlefont = (12, "serif"), 
    guidefont = (10, "serif"), 
    legendfont = (8, "serif"), 
    tickfont = (8, "serif"), 
    # Grids
    minorticks=false,
    minorgrid=false,
    grid=false
    # Image
)

include("FittingUtils.jl")

paramLabels = Dict(
    "K" => ["K", "E", "a", "h", "θ", "K2", "a2"], 
    "J" => ["K", "E", "a", "h", "θ", "α13", "ϵ3", "K2", "a2"]
)

paramLimits = Dict(
    "K" => (0, 100), 
    "E" => (6.4, 7.), 
    "a" => (0, 0.998), 
    "h" => (1, 15), 
    "θ" => (5, 85), 
    "α13" => (-0.4, 10), 
    "ϵ3" => (-0.4, 10), 
    "K2" => (0, 100), 
    "a2" => (0, 10)
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

    stats .- sum(result.stats)

end

function SingleCombo(result, idx1::Int, idx2::Int; paramLimits=paramLimits, paramLabels=paramLabels["J"], num=51)

    label1 = paramLabels[idx1]
    label2 = paramLabels[idx2]

    param1 = range(paramLimits[label1]..., num)
    param2 = range(paramLimits[label2]..., num)

    grid = ChiSquaredGrid(result, idx1, idx2, param1, param2)

    # Finding the standard deviation to plot 1σ,2σ,3σ contours
    stdev = std(grid)
    levels = [stdev, 2*stdev, 3*stdev]

    ctr = contour(param1, param2, grid, levels=[levels[1]], c=[:black], cbar=false, ls=:solid)
    contour!(ctr, param1, param2, grid, levels=[levels[2]], c=[:black], cbar=false, ls=:dash)
    contour!(ctr, param1, param2, grid, levels=[levels[3]], c=[:black], cbar=false, ls=:dot)
    scatter!(ctr, (result.u[idx1], result.u[idx2]), c=:black, msw=0, label=nothing, ticks=nothing)

end

function SingleParam(result, idx; paramLimits=paramLimits, paramLabels=paramLabels["J"], num=51)

    label = paramLabels[idx]
    param = range(paramLimits[label]..., num)

    values = copy(result.u)

    stats = Array{Float64}(undef, num)

    for (i, p) in enumerate(param)
        values[idx] = p
        stats[i] = measure(ChiSquared(), result, values)
    end

    plt = plot(param, stats)
    vline!(plt, result.u[idx])

end

function GetCombos(number)

    println(number)

    combos = Array{Tuple{Int, Int}}(undef, sum(1:number-1))
    k=1

    for i in number:-1:1
        for j in number:-1:i
            if i !== j
                combos[k] = (i,j)
                k+=1
            end
        end
    end

    combos

end

function CornerPlot(result; paramLimits=paramLimits, paramLabels=paramLabels["J"], num=51)

    len = length(result.u)

    combos = GetCombos(len)

    plots = Array{Plots.Plot{Plots.PyPlotBackend}}(undef, len, len)

    for combo in combos
        plt = SingleCombo(result, combo...; paramLimits=paramLimits, paramLabels=paramLabels, num=num)
        plots[combo...] = plt
        plots[len-combo[1]+1, len-combo[2]+1] = plot(ticks=nothing, framestyle=:none)
    end

    for i in 1:len
        plots[i,i] = plot(ticks=nothing, framestyle=:none)
    end

    display(plot(vec(plots)..., layout=grid(len,len), dpi=100, size=cm2px.((12,12))))

    plots

end

@load "output/testfit.bson" johannsenResult

plots = CornerPlot(johannsenResult[1]; num=3)

# plt = SingleCombo(johannsenResult[1], 1, 2; num=5)

close("all")