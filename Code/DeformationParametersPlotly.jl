using Gradus
using Plots
using ProgressBars
using LaTeXStrings

plotly()

cm2in(x) = 0.3937008x
cm2px(x) = Int(round(cm2in(x) * 100))

default(titlefont = (12, "serif"), 
    guidefont = (10, "serif"), 
    legendfont = (8, "serif"), 
    tickfont = (8, "serif"), 
    colorbar_titlefont = (10, "serif"),
    gridalpha=1.,
    minorticks=true,
    dpi=300,
    size=cm2px.((12,9))
)

function ValidityCheck(m)

    if is_naked_singularity(m)
        # Naked singularity
        return -2.
    elseif is_no_isco(m)
        return -1
    else
        return Gradus.isco(m)
    end
end

function is_no_isco(m)
    try 
        Gradus.isco(m)
        return false
    catch
        return true
    end
end

Constraints(a) = -(1+sqrt(1 - a^2))^3

function ParameterRegions(αmax, ϵmax; a=0.998, step=0.1)

    ticks = -10:2:10

    minVal = Constraints(a)

    αs = minVal:step:αmax
    ϵs = minVal:step:ϵmax

    regions = zeros(Float64, (length(αs), length(ϵs)))

    pbar = ProgressBar(total=length(αs)*length(ϵs), printing_delay=0.1)

    Threads.@threads for i in eachindex(αs)
        α = αs[i]
        for (j, ϵ) in enumerate(ϵs)

            update(pbar)
            
            m = JohannsenMetric(M=1., a=a, ϵ3=ϵ, α13=α)
            
            regions[i,j] = ValidityCheck(m)

        end
    end

    hmp = heatmap(
        ϵs,
        αs,
        regions;
        xlabel = "epsilon_3",
        ylabel = "alpha_13",
        clims=(0, maximum(regions)),
        ylims=(minimum(αs), maximum(αs)),
        xlims=(minimum(ϵs), maximum(ϵs)),
        colorbar_title="ISCO (R_g)",
        framestyle=:box,
        title="a = $a",
        minorticks=4,
        minorgrid=true,
        minorgridalpha=0.5,
        aspect_ratio=:equal,
        xticks=ticks[ticks.>minVal],
        yticks=ticks[ticks.>minVal]
    )

    display(hmp)

    return regions, αs, ϵs, hmp
end

function DrawHorizon(p, m)
    rs, θs = event_horizon(m, resolution = 200)
    radius = rs

    x = @. radius * sin(θs)
    y = @. radius * cos(θs)
    plot!(p, x, y, label = "a = $(m.a)")
end
#range(0, 0.998, 5)

for a in range(0, 0.998, 10)
    ParameterRegions(10., 10.; step=0.1, a=a)
end


