include("utils/GeneralUtils.jl")
include("utils/DeformUtils.jl")

using BSON: @save, @load

total=31

as = range(-0.998, 0.998, total)

function BoundPlot(as, pbar)
    for (i, a) in enumerate(as)

        update(pbar)
        regions, _, _, _ = ParameterRegions(10., 10.; a=a, verbose=false)

        @save "outBin/$(round(a; digits=3))" regions

        GC.gc()
        
    end
end

pbar = ProgressBar(total=total)

chunks = Iterators.partition(as, cld(length(as), Threads.nthreads()))
tasks = map(chunks) do chunk
    Threads.@spawn BoundPlot(chunk, pbar)
end
chunk_sums = fetch.(tasks)

for (i, a) in enumerate(as)

    @load "outBin/$(round(a; digits=3))" regions

    minVal = Constraints(a)
    αs = minVal:0.1:10
    ϵs = minVal:0.1:10

    hmp = PlotBounds((regions, αs, ϵs, minVal), a)

    savefig(hmp, "$i.png")

end