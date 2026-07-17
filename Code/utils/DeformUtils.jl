include("Deformations.jl")
include("PlottingDefaults.jl")

using DataFrames, CSV, ProgressBars, Gradus, LaTeXStrings

const DeformationBoundsTable = DataFrame(CSV.File("Code/utils/NewDeformationBounds.csv"))

export FindTopVertices, FindBottomVertices, FindWedgeVertices, FindVertices, FindShape, FindShapes, PlotBounds, CheckBounds, GenerateBoundsCSV, GenerateRegionFiles, PlotError, ParameterRegions, PlotRegion, DeformationSpinPlot, DeformationSpin

using BSON: @save, @load

# =============================================================================
# Boundary Finding
# =============================================================================

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

# =============================================================================
# Boundary Finding
# =============================================================================

function GenerateBoundsCSV(as=range(-0.998, 0.998, 20))

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
    end

    CSV.write("out.csv", df)
end

function GenerateRegionFiles(as=range(-0.998, 0.998, 20))

    for a in as
        regions, _, _ = ParameterRegions(10, 10; a=a, step=0.05)
        @save "outBin/a=$(round(a; digits=3)).bson" regions
    end

end

function PlotBounds(a, df=DeformationBoundsTable; sgns=[1,1,-1,-1,-1,-1], verbose=true)

    plts = ParameterRegions(10, 10; a=a, step=0.4, verbose=verbose)
    plts[1][plts[1] .== -3] .+= 1
    hmp = PlotRegion(plts...; title="a=$(round(a; digits=3))", c=:grays)

    lines = GetLines(a, df)

    tri1Int = (lines[2].c - lines[1].c) / (lines[1].m - lines[2].m)
    tri2Int = (lines[4].c - lines[3].c) / (lines[3].m - lines[4].m)

    x = minimum(plts[3]):0.05:maximum(plts[3])
    
    ys = Array{Array{Float64}}(undef, 4)

    fill = 10*ones(length(x))
    fillranges = (-fill, -fill, fill, -fill)
    fillstyles = (:..., :..., :..., :...)

    ys[1] = LineEquation(x, lines[5])
    ys[2] = LineEquation(x, lines[6])
    ys[3] = append!(LineEquation(x[x.<=tri1Int], lines[1]), LineEquation(x[x.>=tri1Int], lines[2]))
    ys[4] = append!(LineEquation(x[x.<=tri2Int], lines[3]), LineEquation(x[x.>=tri2Int], lines[4]))
    
    colours = [:gray, :gray, :black, :black]

    for (i, y) in enumerate(ys)
        plot!(
            x, y, 
            label=nothing, 
            c=colours[i], 
            fillrange=fillranges[i],
            fillstyle=fillstyles[i],
            cbar=false
        )
    end

    return hmp
end

function CheckBounds(N; df=DeformationBoundsTable, as=-0.998:0.0001:0.998)

    out = nothing

    pbar = ProgressBar(total=N)

    Threads.@threads for a in rand(as, N)

        update(pbar)
        
        α13 = rand(Constraints(a):0.0001:10)
        ϵ3 = rand(Constraints(a):0.0001:10)

        m = JohannsenMetric(a=a, α13 = α13, ϵ3 = ϵ3)

        if QuickIsValid(ϵ3, α13, a, df) & !IsValid(ϵ3, α13, a)
            println("Combination ϵ3=$ϵ3, α13=$α13, a=$a failed")

            m = JohannsenMetric(ϵ3=ϵ3, α13=α13, a=a)
            println(ValidityCheckISCO(m))
            println()

            out = (ϵ3, α13, a)
            break
        end

    end
    
    GC.gc()

    return out

end

function PlotError(ϵ3, α13, a; df=DeformationBoundsTable)
    hmp = PlotBounds(a, df)
    scatter!(hmp, [ϵ3], [α13]; c=:white, label=nothing, msc=:red)
    display(hmp)
end

function ParameterRegions(αmax, ϵmax; a=0.998, step=0.1, verbose=true)
    """
    Checking the validity for a grid of points in the parameter space for a given spin value and producing a heatmap
    """
    minVal = Constraints(a)

    αs = minVal:step:αmax
    ϵs = minVal:step:ϵmax

    regions = zeros(Float64, (length(αs), length(ϵs)))
    
    if verbose
        pbar = ProgressBar(total=length(αs)*length(ϵs), printing_delay=0.1)
    end

    for i in eachindex(αs)
        α = αs[i]
        for (j, ϵ) in enumerate(ϵs)

            if verbose
                update(pbar)
            end
            
            m = JohannsenMetric(M=1., a=a, ϵ3=ϵ, α13=α)
            
            regions[i,j] = ValidityCheckISCO(m)

        end
    end

    return regions, αs, ϵs, minVal
end

function PlotRegion(regions, αs, ϵs, minVal; ticks = -10:2:10, title="", contour=false, c=:default)
    """
    Plotting the output from ParameterRegions
    """

    hmp = heatmap(
        ϵs,
        αs,
        regions;
        xlabel = L"\epsilon_3",
        ylabel = L"\alpha_{13}",
        clims=(minimum(regions), 0),
        ylims=(minimum(αs), maximum(αs)),
        xlims=(minimum(ϵs), maximum(ϵs)),
        colorbar_title=L"ISCO $(R_g)$",
        title=title,
        minorticks=4,
        minorgrid=true,
        minorgridalpha=0.5,
        aspect_ratio=:equal,
        xticks=ticks[ticks.>minVal],
        yticks=ticks[ticks.>minVal],
        c=c
    )

end

function DeformationSpinPlot(x, y, regions; xlabel=L"a", ylabel=L"\epsilon_3")
    """
    Plotting the output from DeformationSpin
    """

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

function DeformationSpin(;xs = -0.998:0.005:0.998, ys = -10:0.1:10, param="α13", plt=true)
    """
    Checking the 2D parameter space of a deformation parameter and black hole spin
    for validity and plotting.
    This produces the same result as Figure 6 in Johannsen (2013)
    """

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
            
            regions[i,j] = ValidityCheckISCO(m)

        end
    end

    regions[regions.<0] .= NaN

    if plt
        Figure6(xs, ys, regions; ylabel=param)
    end

    xs, ys, regions

end