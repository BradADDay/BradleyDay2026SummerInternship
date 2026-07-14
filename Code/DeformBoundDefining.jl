using BSON: @save, @load
using DataFrames
using Flatten: flatten
using CSV
using Interpolations

include("utils/UTILS.jl")

function Gradient(coord1, coord2)
    diff = coord1 .- coord2
    return diff[2]/diff[1]
end

function Intercept(grad, coord)
    return coord[2] - grad * coord[1]
end

struct Line{T<:Real}
    "Gradient"
    m::T
    "Intercept"
    c::T
end

function Line(coord1::AbstractArray, coord2::AbstractArray)
    m = Gradient(coord1, coord2)
    c = Intercept(m, coord1)
    Line(m, c)
end

function LineEquation(x, line::Line)
    line.m * x .+ line.c
end

function FindTopVertices(regions, ϵs, αs)

    len=length(regions[1,:])

    vertices = Array{Any}(undef, 3, 2)

    if regions[end, 1] == -3
        array = regions[:,1]
        vertices[1,:] = [1, len-findfirst(!=(-3), reverse(array))]
    else
        array = regions[end,:]
        vertices[1,:] = [findfirst(==(-3), array), len]
    end

    found = true
    i = 0

    while found
        found = false
        idx = findfirst(==(-3), regions[end-i,:])
        if !isnothing(idx)
            vertices[2,:] = [idx, len-i]
            found = true
        end

        i+=1
    end

    begin
        array = regions[end,:]
        x = findlast(==(-3), array)
        vertices[3,1] = x
        if x != len
            vertices[3,2] = len
        else
            y = len - findlast(==(-3), reverse(regions[:,end]))
            vertices[3,2] = y
        end
    end
    vertices[:,1] .= ϵs[vertices[:,1]]
    vertices[:,2] .= αs[vertices[:,2]]
    
    vertices

end

function FindBottomVertices(regions, ϵs, αs)

    len=length(regions[1,:])

    vertices = Array{Any}(undef, 3, 2)

    if regions[1, 1] == -3
        array = regions[:,1]
        vertices[1,:] = [1, findfirst(!=(-3), array)]
    else
        array = regions[1,:]
        vertices[1,:] = [findfirst(==(-3), array), 1]
    end

    found = true
    i = 1

    while found
        found = false
        idx = findfirst(==(-3), regions[i,:])

        if !isnothing(idx)
            vertices[2,:] = [idx, i]
            found = true
        end

        i+=1
    end

    begin
        array = regions[1,:]
        vertices[3,:] = [findlast(==(-3), array), 1]
    end

    vertices[:,1] .= ϵs[vertices[:,1]]
    vertices[:,2] .= αs[vertices[:,2]]
    
    vertices

end

function FindWedgeVertices(regions, ϵs, αs)

    if regions[1, end] == -1

        len=length(regions[1,:])

        vertices = Array{Any}(undef, 2, 2)

        begin
            array = regions[:,end]
            vertices[1,:] = [len, findlast(==(-1), array)]
        end

        begin
            array = regions[1, :]
            vertices[2,:] = [findfirst(==(-1), array), 1]
        end

        vertices[:,1] .= ϵs[vertices[:,1]]
        vertices[:,2] .= αs[vertices[:,2]]
        
        vertices
    end
end

function FindVertices(regions, ϵs, αs; verbose=false)

    top = FindTopVertices(regions, ϵs, αs)

    bottom = FindBottomVertices(regions, ϵs, αs)

    wedge = FindWedgeVertices(regions, ϵs, αs)

    if verbose
        if isnothing(top)
            println("WARNING: No top triangle found.")
        end
        if isnothing(bottom)
            println("WARNING: No bottom triangle found.")
        end
        if isnothing(wedge)
            println("WARNING: No wedge found.")
        end
    end

    top, bottom, wedge
end

function FindShape(vertices, numLines)

    if !isnothing(vertices)
        out = Array{Any}(undef, numLines)

        for i in 1:numLines
            m = Gradient(vertices[i,:], vertices[i+1,:])
            out[i] = m, Intercept(m, vertices[i,:])
        end

        return out
    else 
        return [nothing, nothing]
    end
end

function FindShapes(vertices)

    top1, top2 = FindShape(vertices[1], 2)
    bottom1, bottom2 = FindShape(vertices[2], 2)
    wedge = FindShape(vertices[3], 1)[1]

    top1, top2, bottom1, bottom2, wedge
    
end

#= as = range(-0.998, 0.998, 20)

# for a in as
#     regions, _, _ = ParameterRegions(10, 10; a=a, step=0.05)
#     @save "outBin/a=$(round(a; digits=3)).bson" regions
# end

a = as[3]

cols = 

df = DataFrame([[],[],[],[],[],[],[],[],[],[],[],[],[]], ["a", "m1", "c1", "m2", "c2", "m3", "c3", "m4", "c4", "m5", "c5", "m6", "c6"])

for a in as
    @load "outBin/a=$(round(a; digits=3)).bson" regions
        
    minVal = Constraints(a)

    αs = minVal:0.05:10
    ϵs = minVal:0.05:10

    vertices = FindVertices(regions, ϵs, αs)
    shapes = FindShapes(vertices)

    line = Dict{String,Any}()

    line["a"] = a

    for i in eachindex(shapes)
        if !isnothing(shapes[i])
            line["m$i"] = shapes[i][1]
            line["c$i"] = shapes[i][2]
        else
            line["m$i"] = NaN
            line["c$i"] = NaN
        end
    end

    line["m6"] = NaN
    line["c6"] = NaN

    push!(df, line)

    hmp = PlotRegion(regions, αs, ϵs, minVal; title=a)

    for shape in shapes

        if !isnothing(shape)
            plot!(ϵs, LineEquation(ϵs, Line(shape...)); c=:grey)
        end

    end
    display(hmp)
end =#

# CSV.write("test.csv", df)

df = DataFrame(CSV.File("test.csv"))

function GetLines(a, df)

    df = copy(df)

    as = df.a
    data = select!(copy(df), Not(:a))

    len = length(names(data))

    out = Array{Float64}(undef, len)

    for (i,name) in enumerate(String.(names(data)))
        interp = linear_interpolation(as, data[!, name])
        out[i] = interp(a)
    end

    returns = Array{Any}(undef, Int(round(len/2)))

    for i in 1:Int(round(len/2))
        if !isnan(out[2*i-1])
            returns[i] = Line(out[2*i-1], out[2*i]) 
        else
            returns[i] = nothing
        end
    end

    returns
end

function PlotBounds(a, df; padding=0.1, sgns=[1,1,-1,-1,-1,-1])
    plts = ParameterRegions(10, 10; a=a, step=0.5)
    hmp = PlotRegion(plts...; title="$a")

    lines = GetLines(a, df)

    for (i, line) in enumerate(lines)
        if !isnothing(line)
            plot!(plts[3], LineEquation(plts[3], line) .- padding*sgns[i], label="$i")
        end
    end

    display(hmp)
    return hmp
end

N = 100_000

pbar = ProgressBar(total=N)

Threads.@threads for a in rand(-0.998:0.001:0.998, N)
    update(pbar)
    α13 = rand(Constraints(a):0.001:10)
    ϵ3 = rand(Constraints(a):0.001:10)

    m = JohannsenMetric(a=a, α13 = α13, ϵ3 = ϵ3)

    if QuickIsValid(ϵ3, α13, a, df) & !IsValid(ϵ3, α13, a)
        println("Combination ϵ3=$ϵ3, α13=$α13, a=$a failed")
        println()
        break
    end

end

close("all")
GC.gc()

##

ϵ3=-2.628753736454565
α13=-1.306753736454565
a=-0.92

hmp = PlotBounds(a, df)
scatter!(hmp, [ϵ3], [α13]; c=:white, label=nothing, msc=:red)
display(hmp)