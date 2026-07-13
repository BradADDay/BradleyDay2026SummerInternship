using BSON: @save, @load

include("utils/DeformUtils.jl")
include("utils/PlottingDefaults.jl")

function FindTopVertices(regions, αs, ϵs)

    len=length(regions[1,:])

    vertices = Array{Any}(undef, 3, 2)

    println("Finding vertex 1")
    if regions[end, 1] == -3
        println("Searching Column")
        array = regions[:,1]
        vertices[1,:] = [1, len-findfirst(!=(-3), reverse(array))]
    else
        println("Searching Row")
        array = regions[end,:]
        vertices[1,:] = [findfirst(==(-3), array), len]
    end

    println("Finding vertex 2")
    begin
        array = regions[end,:]
        x = findlast(==(-3), array)
        vertices[2,1] = x
        if x != len
            vertices[2,2] = len
        else
            y = len - findlast(==(-3), reverse(regions[:,end]))
            vertices[2,2] = y
        end
    end

    found = true
    i = 0

    println("Finding vertex 3")
    while found
        found = false
        idx = findfirst(==(-3), regions[end-i,:])
        if !isnothing(idx)
            vertices[3,:] = [idx, len-i]
            found = true
        end

        i+=1
    end

    vertices[:,1] .= ϵs[vertices[:,1]]
    vertices[:,2] .= αs[vertices[:,2]]
    
    vertices

end

function FindBottomVertices(regions, αs, ϵs)

    len=length(regions[1,:])

    vertices = Array{Any}(undef, 3, 2)

    println("Finding vertex 1")
    if regions[1, 1] == -3
        println("Searching Column")
        array = regions[:,1]
        vertices[1,:] = [1, findfirst(!=(-3), array)]
    else
        println("Searching Row")
        array = regions[1,:]
        vertices[1,:] = [findfirst(==(-3), array), 1]
    end

    println("Finding vertex 2")
    begin
        array = regions[1,:]
        vertices[2,:] = [findlast(==(-3), array), 1]
    end

    found = true
    i = 1

    println("Finding vertex 3")
    while found
        found = false
        idx = findfirst(==(-3), regions[i,:])

        if !isnothing(idx)
            vertices[3,:] = [idx, i]
            found = true
        end

        i+=1
    end

    vertices[:,1] .= ϵs[vertices[:,1]]
    vertices[:,2] .= αs[vertices[:,2]]
    
    vertices

end

function FindWedgeVertices(regions, αs, ϵs)

    if regions[1, end] == -1

        len=length(regions[1,:])

        vertices = Array{Any}(undef, 2, 2)

        println("Finding vertex 1")
        begin
            array = regions[:,end]
            vertices[1,:] = [len, findlast(==(-1), array)]
        end

        println("Finding vertex 2")
        begin
            array = regions[1, :]
            vertices[2,:] = [findfirst(==(-1), array), 1]
        end

        vertices[:,1] .= ϵs[vertices[:,1]]
        vertices[:,2] .= αs[vertices[:,2]]
        
        vertices
    end
end

as = range(-0.998, 0.998, 20)

# for a in as
#     regions, _, _ = ParameterRegions(10, 10; a=a, step=0.05)
#     @save "outBin/a=$(round(a; digits=3)).bson" regions
# end

a = as[4]

@load "outBin/a=$(round(a; digits=3)).bson" regions
    
minVal = Constraints(a)

αs = minVal:0.05:10
ϵs = minVal:0.05:10

hmp = PlotRegion(regions, αs, ϵs, minVal; title=a)
vertices = FindTopVertices(regions, αs, ϵs)
scatter!(vertices[:,1], vertices[:,2]; label=nothing, markersize=3, c=:white, msw=0.2)
vertices = FindBottomVertices(regions, αs, ϵs)
scatter!(vertices[:,1], vertices[:,2]; label=nothing, markersize=3, c=:white, msw=0.2)
vertices = FindWedgeVertices(regions, αs, ϵs)
scatter!(vertices[:,1], vertices[:,2]; label=nothing, markersize=3, c=:white, msw=0.2)
