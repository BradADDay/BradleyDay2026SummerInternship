include("utils/GeneralUtils.jl")
include("utils/DeformUtils.jl")
include("utils/PlottingDefaults.jl")

using .BoundsDefinition

N = 5000

out = nothing
i=0

while true
    
    out = BoundsDefinition.CheckBounds(N)

    if !isnothing(out)
        PlotError(out...)

        close("all")
        GC.gc()
        CompleteSound()
        break
    end

    i+=N

    if i > 1E6
        print("No bad values found")
        break
    end

end

