"""
File containing a selection of utility functions for investigating the deformation parameter space
"""

include("Deformations.jl")
include("PlottingDefaults.jl")

using ProgressBars, LaTeXStrings, JSON3
using BSON: @save, @load

# =============================================================================
# Boundary Finding
# =============================================================================

function FindTopVertices(regions, ϵs, αs)
    """
    Find the vertices of the top triangle in the parameter space plot
    """

    # The size of the grid
    len=length(regions[1,:])

    vertices = Array{Any}(undef, 3, 2)

    # Finding the first corner (top left)
    # -----------------------------------
    # This finds the first invalid pixel from the top left. 
    # If the top left pixel is invalid then it searches the first column 
    #    to find the last invalid pixel from the top

    if regions[end, 1] == -3
        array = regions[:,1]
        vertices[1,:] = [1, len-findfirst(!=(-3), reverse(array))]
    else
        array = regions[end,:]
        vertices[1,:] = [findfirst(==(-3), array), len]
    end

    # Finding the second corner (bottom middle)
    # -----------------------------------------
    # This scans row by row from the top left, seeing if an invalid pixel is present, 
    #    if there is no invalid pixel in a given row then the loop is broken and 
    #    the most recent invalid pixel is returned

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

    # Finding the third corner (top right)
    # ------------------------------------
    # This finds the last invalid pixel in the top row
    # If the top right pixel is invalid, it searches the last column to find where 
    #    the invalid region intersects the leftmost edge of the image

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

    # Converting pixel coordinates to parameter space coordinates
    vertices[:,1] .= ϵs[vertices[:,1]]
    vertices[:,2] .= αs[vertices[:,2]]
    
    vertices

end

function FindBottomVertices(regions, ϵs, αs)
    """
    Find the vertices of the bottom triangle in the parameter space plot
    """

    vertices = Array{Any}(undef, 3, 2)

    # Finding the first corner (top left)
    # -----------------------------------
    # This finds the first invalid pixel from the bottom left. 
    # If the bottom left pixel is invalid then it searches the first column 
    #    to find the last invalid pixel from the bottom.

    if regions[1, 1] == -3
        array = regions[:,1]
        vertices[1,:] = [1, findfirst(!=(-3), array)]
    else
        array = regions[1,:]
        vertices[1,:] = [findfirst(==(-3), array), 1]
    end

    # Finding the second corner (top middle)
    # -----------------------------------------
    # This scans row by row from the bottom left, seeing if an invalid pixel is present, 
    #    if there is no invalid pixel in a given row then the loop is broken and 
    #    the most recent invalid pixel is returned

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

    # Finding the third corner (bottom right)
    # ------------------------------------
    # This finds the last invalid pixel in the bottom row

    begin
        array = regions[1,:]
        vertices[3,:] = [findlast(==(-3), array), 1]
    end

    # Converting pixel coordinates to parameter space coordinates
    vertices[:,1] .= ϵs[vertices[:,1]]
    vertices[:,2] .= αs[vertices[:,2]]
    
    vertices

end

function FindWedgeVertices(regions, ϵs, αs)
    """
    Find the vertices of the bottom right wedge shape
    """

    # Checking if the bottom right pixel is invalid
    if regions[1, end] == -1

        # Size of the parameter grid
        len=length(regions[1,:])

        vertices = Array{Any}(undef, 2, 2)

        # Finding the last invalid pixel in the rightmost column
        begin
            array = regions[:,end]
            vertices[1,:] = [len, findlast(==(-1), array)]
        end

        # Finding the first invalid pixel in the bottom row
        begin
            array = regions[1, :]
            vertices[2,:] = [findfirst(==(-1), array), 1]
        end

        
        # Converting pixel coordinates to parameter space coordinates
        vertices[:,1] .= ϵs[vertices[:,1]]
        vertices[:,2] .= αs[vertices[:,2]]
        
        vertices
    end
end

function FindVertices(regions, ϵs, αs; verbose=false)
    """
    Find all vertices for invalid regions
    """

    top = FindTopVertices(regions, ϵs, αs)

    bottom = FindBottomVertices(regions, ϵs, αs)

    wedge = FindWedgeVertices(regions, ϵs, αs)

    # Checking vertices were found for all regions
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
    """
    Generate a selection of pixel grids depicting the parameter space
    """

    for a in as
        regions, _, _ = ParameterRegions(10, 10; a=a, step=0.05)
        @save "outBin/a=$(round(a; digits=3)).bson" regions
    end

end

# =============================================================================
# Boundary Validation
# =============================================================================

function PlotBounds(a, df=DeformationBoundsTable; sgns=[1,1,-1,-1,-1,-1], step=0.4, pbarLabel=nothing, verbose=true)
    """
    Plot the parameter space showing the invalid regions, and plot the bounds defined within the table dataframe
    """

    # Plotting the parameter space regions
    plts = ParameterRegions(10, 10; a=a, step=step, verbose=verbose)
    hmp = PlotRegion(plts...; title="a=$(round(a; digits=3))", c=:grays)

    # Getting the boundary lines for the input spin
    lines = GetLines(a, df)

    # Getting the intersections for the top and bottom lines
    tri1Int = Intersection(lines[1], lines[2])
    tri2Int = Intersection(lines[3], lines[4])

    # Getting an array for the x positions
    x = minimum(plts[3]):0.05:maximum(plts[3])
    
    ys = Array{Array{Float64}}(undef, 4)

    # Defining the lines as coordinates for plotting
    ys[1] = LineEquation(x, lines[5])
    ys[2] = LineEquation(x, lines[6])
    ys[3] = append!(LineEquation(x[x.<=tri1Int], lines[1]), LineEquation(x[x.>=tri1Int], lines[2]))
    ys[4] = append!(LineEquation(x[x.<=tri2Int], lines[3]), LineEquation(x[x.>=tri2Int], lines[4]))

    # Defining fill ranges for the Invalid regions
    fill = 10*ones(length(x))
    fillranges = (-fill, -fill, fill, -fill)
    fillstyles = (:..., :..., :..., :...)
    
    # Defining the colours for the lines
    colours = [:gray, :gray, :black, :black]

    # Plotting
    for (i, y) in enumerate(ys)
        plot!(hmp,
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
    """
    Test the defined boundary lines for N random points in the parameter space
    """

    out = nothing

    pbar = ProgressBar(total=N)

    # Looping through a random selection of spin values
    Threads.@threads for a in rand(as, N)

        update(pbar)
        
        # Getting the random deformation parameters
        α13 = rand(Constraints(a):0.0001:10)
        ϵ3 = rand(Constraints(a):0.0001:10)

        if QuickIsValid(ϵ3, α13, a, df) & !IsValid(ϵ3, α13, a)

            println("Combination ϵ3=$ϵ3, α13=$α13, a=$a failed")

            out = (ϵ3, α13, a)

            PlotError(ϵ3, α13, a)

            break
        end

    end
    
    GC.gc()

    return out

end

function PlotError(ϵ3, α13, a; df=DeformationBoundsTable)
    """
    Plot the failed combination of parameters in parameter space
    """

    # Plotting the bounds for inspection
    hmp = PlotBounds(a, df)

    # Plotting the point that failed the checks
    scatter!(hmp, [ϵ3], [α13]; c=:white, label=nothing, msc=:red)
    display(hmp)
end

function ParameterRegions(αmax, ϵmax; a=0.998, step=0.1, pbarLabel=nothing, verbose=true)
    """
    Checking the validity for a grid of points in the parameter space for a given spin value and producing a heatmap
    """
    # Get the minimum valid value using Eqs 75, 77 of Johannsen 2013
    minVal = Constraints(a)

    # Define the parameter space
    αs = minVal:step:αmax
    ϵs = minVal:step:ϵmax

    # Initialise the regions matrix
    regions = zeros(Float64, (length(αs), length(ϵs)))
    
    # Starting a progress bar
    if verbose
        pbar = ProgressBar(total=length(αs)*length(ϵs), printing_delay=0.1)
        set_description(pbar, pbarLabel)
    end

    # Checking the validity for each point on the grid
    for (i, α) in enumerate(αs)
        for (j, ϵ) in enumerate(ϵs)
            
            m = JohannsenMetric(M=1., a=a, ϵ3=ϵ, α13=α)
            
            regions[i,j] = ValidityCheckISCO(m)

            if verbose
                update(pbar)
            end

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

function FixedParamsJSON(
    outname; 
    as=range(-0.998, 0.998, 13),
    α13s=range(-8., 10., 10),
    ϵ3s=range(-8., 10., 10)
    )
    """
    Generate a JSON file containing the corrected parameter space coordinates
    """

    dict = Dict()

    pbar = ProgressBar(total=length(as)*length(α13s)*length(ϵ3s))
    set_description(pbar, "Generating JSON: ")

    for a in as
        dict["$(round(a; digits=3))"] = Dict()
        for α13 in α13s
            oα13 = α13
            for ϵ3 in ϵ3s

                α13 = oα13
                name = "$ϵ3, $α13"

                update(pbar)

                # Finding the corrected coordinates
                if !QuickIsValid(ϵ3, α13, a)
                    ϵ3, α13 = FindNearestSafePoint([ϵ3, α13], a)
                end

                # Storing to dict
                dict["$(round(a; digits=3))"][name] = ϵ3, α13

            end
        end
    end

    # Saving to JSON
    open("$outname.json", "w") do f
        JSON3.write(f, dict)
        println(f)
    end

end

function CheckDeformCorrection(;
    as = range(-0.998, 0.998, 10),
    ϵs = range(-8, 10, 50),
    αs = range(-8, 10, 50)
    )
    """
    Generate a grid of parameter space coordinates and plot where they are 
        corrected to by FindNearestSafePoint
    """

    for a in as

        # Plotting the boundaries
        plt = PlotBounds(a)

        # Checking if a point on the grid is valid, correcting if not
        for ϵ in ϵs, α in αs
            if !QuickIsValid(ϵ, α, a)
                ϵ, α = FindNearestSafePoint([ϵ, α], a)
                scatter!(plt, [ϵ], [α]; c=:red, msw=0, label=nothing, markersize=2)
            end
        end
        
        display(plt)
        
    end
end

# =================================================================================
# Figure 6 of Johannsen 2013
# =================================================================================

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

