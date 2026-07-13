include("utils/DeformUtils.jl")

as = range(-0.998, 0.998, 20)

for a in as
    regions, _, _ = ParameterRegions(10, 10; a=a, step=0.05)
    @save "outBin/a=$(round(a; digits=3)).bson" regions
    GC.gc()
end

for a in as
    @load "outBin/a=$(round(a; digits=3)).bson" regions
    
    minVal = Constraints(a)

    αs = minVal:0.2:10
    ϵs = minVal:0.2:10

    PlotRegion(regions, αs, ϵs, minVal)

end