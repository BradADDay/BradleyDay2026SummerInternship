include("utils/GeneralUtils.jl")
include("utils/DeformUtils.jl")

total=31

as = range(-0.998, 0.998, total)

function BoundPlot(as, pbar)
    for (i, a) in enumerate(as)

        update(pbar)
        plt = PlotBounds(a; step = 0.2, verbose=false)

        savefig(plt, "outBin/$i_New.png")

        close("all")
        GC.gc()
        
    end
end

pbar = ProgressBar(total=total)

chunks = Iterators.partition(as, cld(length(as), Threads.nthreads()))
tasks = map(chunks) do chunk
    Threads.@spawn BoundPlot(chunk, pbar)
end
chunk_sums = fetch.(tasks)