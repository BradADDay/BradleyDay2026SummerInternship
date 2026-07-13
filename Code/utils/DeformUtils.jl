using Gradus
using ProgressBars
using LaTeXStrings

include("PlottingDefaults.jl")

function ValidityCheck(m)

    if is_naked_singularity(m)
        # Naked singularity
        return -3
    elseif ((m.α13 < Constraints(m.a)) | (m.ϵ3 < Constraints(m.a)))
        # Abnormal exterior region
        return -2
    elseif is_no_isco(m)
        # No ISCO
        return -1
    else
        # No Abnormalities
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

    return regions, αs, ϵs
end

function DrawHorizon(p, m)
    rs, θs = event_horizon(m, resolution = 200)
    radius = rs

    x = @. radius * sin(θs)
    y = @. radius * cos(θs)
    plot!(p, x, y, label = "a = $(m.a)")
end

function PlotRegion(regions, αs, ϵs, minVal; ticks = -10:2:10, title="")

    hmp = heatmap(
        ϵs,
        αs,
        regions;
        xlabel = L"\epsilon_3",
        ylabel = L"\alpha_{13}",
        clims=(-3, 0),
        ylims=(minimum(αs), maximum(αs)),
        xlims=(minimum(ϵs), maximum(ϵs)),
        colorbar_title=L"ISCO $(R_g)$",
        title=title,
        minorticks=4,
        minorgrid=true,
        minorgridalpha=0.5,
        aspect_ratio=:equal,
        xticks=ticks[ticks.>minVal],
        yticks=ticks[ticks.>minVal]
    )

    contour!(ϵs, αs, regions, levels=[-3], c=[:red], lw=1)

    return hmp
end

function Figure6(x, y, regions; xlabel=L"a", ylabel=L"\epsilon_3")

    redbounds = Constraints.(x)
    redbounds[redbounds.<-5] .= -5

    grnup = Constraints.(x)
    grnup[grnup.<-5] .= -5
    grndown = Constraints.(x)

    plot(x, grnup, label=nothing, fillstyle=:xxxx, c=:green, fillrange=grndown)
    plot!(x, zeros(length(x)).+10, label=nothing, c=:red, fillalpha=0.5, fillrange=redbounds)
    heatmap!(
        x,
        y,
        regions;
        xlabel = xlabel,
        ylabel = ylabel,
        clims=(0, maximum(regions)),
        xlims=(minimum(x), maximum(x)),
        ylims=(minimum(y), maximum(y)),
        colorbar_title=L"ISCO $(R_g)$",
        title="",
        c=:grays
    )

    plot!(x, Constraints.(x), label=nothing, fillstyle=://, c=:black, fillrange=zeros(length(x)).-10)
    contour!(x, y, regions, levels=1:12, c=[:red], lw=1, clabels=true)
end

function DeformSpin(;xs = -0.998:0.005:0.998, ys = -10:0.1:10, param="α13")

    regions = zeros(Float64, (length(ys), length(xs)))

    pbar = ProgressBar(total=length(xs)*length(ys), printing_delay=0.1)

    Threads.@threads for i in eachindex(ys)
        y = ys[i]
        for (j, x) in enumerate(xs)

            update(pbar)
            
            if param == "α13"
                m = JohannsenMetric(M=1., a=x, α13=y)
            elseif param == "ϵ3"
                m = JohannsenMetric(M=1., a=x, ϵ3=y)
            end
            
            regions[i,j] = ValidityCheck(m)

        end
    end

    regions[regions.<0] .= NaN

    Figure6(xs, ys, regions; ylabel=param)

end