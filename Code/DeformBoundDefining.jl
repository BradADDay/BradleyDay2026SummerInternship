include("utils/GeneralUtils.jl")
include("utils/DeformUtils.jl")
include("utils/PlottingDefaults.jl")

df = DataFrame(CSV.File("Code/utils/NewDeformationBounds.csv"))

as = range(-0.998, 0.998, 500)

#= open("Code/utils/NewDeformationBounds.csv", "w") do io
    write(io, "a,m1,c1,m2,c2,m3,c3,m4,c4,m5,c5,m6,c6\n")
end

for a in as
    GetLineParams(a, df)

    open("Code/utils/NewDeformationBounds.csv", "a") do io
        write(io, "$(round(a; digits=3)), $(round.(GetLineParams(a, df)[1]; digits=3))\n")
    end
end =#

pbar = ProgressBar(total=10)

for a in range(0.12, 0.16, 10)[1]

    update(pbar)

    plot = PlotBounds(a, df; verbose=false)
    display(plot)

    savefig(plot, "outBin/$a.png")

    GC.gc()
    close("all")

end


ϵ3, α13, a = CheckBounds(100000; as=-0.96:0.0001:-0.92)

df = DataFrame(CSV.File("Code/utils/NewDeformationBounds.csv"))

ϵ3=1.239810678917241
α13=-0.712189321082759
a=-0.9443

PlotError(ϵ3, α13, a; df)
##

include("utils/GeneralUtils.jl")
include("utils/DeformUtils.jl")
include("utils/PlottingDefaults.jl")

as = range(-0.998, 0.998, 10)
ϵs = range(-8, 10, 50)
αs = range(-8, 10, 50)

for a in as

    plt = PlotBounds(a)

    for ϵ in ϵs, α in αs
        if !QuickIsValid(ϵ, α, a)
            ϵ, α = FindNearestSafePoint([ϵ, α], a)
            scatter!(plt, [ϵ], [α]; c=:red, msw=0, label=nothing, markersize=2)
        end
    end
    xlims!(-9, 15)
    ylims!(-9, 11)
    display(plt)
end

##

using JSON3

include("utils/GeneralUtils.jl")
include("utils/DeformUtils.jl")
include("utils/PlottingDefaults.jl")

as   = range(-0.998, 0.998, 13)
α13s = range(-8., 10., 10)
ϵ3s  = range(-8., 10., 10)

dict = Dict()

pbar = ProgressBar(total=length(as)*length(α13s)*length(ϵ3s))

for a in as
    dict["$(round(a; digits=3))"] = Dict()
    for α13 in α13s
        oα13 = α13
        for ϵ3 in ϵ3s

            α13 = oα13
            name = "$ϵ3, $α13"
            println(a, " ", α13, " ", ϵ3)

            update(pbar)

            if !QuickIsValid(ϵ3, α13, a)
                ϵ3, α13 = FindNearestSafePoint([ϵ3, α13], a)
            end
            println(a, " ", α13, " ", ϵ3)

            dict["$(round(a; digits=3))"][name] = ϵ3, α13

        end
    end
end

open("my.json", "w") do f
    JSON3.write(f, dict)
    println(f)
end

#0.998 4.0
#-2.0 6.0

plt=PlotBounds(0.998)

FindNearestSafePoint([6., -2.], 0.998)