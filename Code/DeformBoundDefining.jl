include("utils/GeneralUtils.jl")
include("utils/DeformUtils.jl")
include("utils/PlottingDefaults.jl")
##
df = DataFrame(CSV.File("Code/utils/NewDeformationBounds.csv"))

open("Code/utils/NewDeformationBounds.csv", "w") do io
    write(io, "a,m1,c1,m2,c2,m3,c3,m4,c4,m5,c5,m6,c6\n")
end

for a in range(-0.998, 0.998, 500)
    GetLineParams(a, df)

    open("Code/utils/NewDeformationBounds.csv", "a") do io
        write(io, "$(round(a; digits=3)), $(round.(GetLineParams(a, df)[1]; digits=3))\n")
    end
end

##

df = DataFrame(CSV.File("Code/utils/NewDeformationBounds.csv"))

for a in range(-0.998, 0.998, 500)[0:50]
    plot = PlotBounds(a, df)

    savefig(plot, "outBin/$a.png")

    GC.gc()
    close("all")
end
