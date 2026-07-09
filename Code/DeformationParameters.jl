using Gradus
using Plots
using ProgressBars
using LaTeXStrings

pyplot()

cm2in(x) = 0.3937008x
cm2px(x) = Int(round(cm2in(x) * 100))

default(titlefont = (12, "serif"), 
    guidefont = (10, "serif"), 
    legendfont = (8, "serif"), 
    tickfont = (8, "serif"), 
    gridalpha=1.,
    minorticks=true,
    dpi=300,
    size=cm2px.((12,9))
)

function ValidityCheck(m)

    if is_naked_singularity(m)
        return NaN
    else
        try
            ISCO = Gradus.isco(m)
            return ISCO
        catch
            return -1.
        end
    end

end

Constraints(a) = -(1+sqrt(1 - a^2))^3

function ParameterRegions(αmax, ϵmax; a=0.998, step=0.1)

    αs = Constraints(a):step:αmax
    ϵs = Constraints(a):step:ϵmax

    regions = zeros(Float64, (length(αs), length(ϵs)))

    pbar = ProgressBar(total=length(αs)*length(ϵs))

    Threads.@threads for i in eachindex(αs)
        α = αs[i]
        for (j, ϵ) in enumerate(ϵs)

            update(pbar)
            
            m = JohannsenMetric(M=1., a=a, ϵ3=ϵ, α13=α)
            
            regions[i,j] = ValidityCheck(m)

        end
    end

    heatmap(
        αs,
        ϵs,
        regions;
        ylabel = L"\epsilon_3",
        xlabel = L"\alpha_{13}",
        clims=(-1, maximum(regions)),
        xlims=(minimum(αs), maximum(αs)),
        ylims=(minimum(ϵs), maximum(ϵs)),
        colorbar_title=L"ISCO $(R_g)$",
        framestyle=:box,
        title="Valid Deformation Parameter Combinations",
        minorticks=4,
        minorgrid=true,
        minorgridalpha=0.5,
        aspect_ratio=:equal
    )
end

ParameterRegions(10., 10.; step=0.02)
