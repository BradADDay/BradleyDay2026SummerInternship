using Gradus
using Plots
using ColorSchemes
using SpectralFitting

# =======================================================================
# Functions
# =======================================================================

function generateConfig(variables, setupDict)
    """
    Generate a string denoting the configuration of a given data run.
    """
    config = ""

    # Looping through the variables
    for i in variables
        config *= "$i = $(round(setupDict[i]; digits= 2))"
        if i != variables[end]
            config *= ", "
        end
    end
    return config
end

function ComputeLineProfile(m, x; height, bins = range(0.0, 1.5, 180), minrₑ=-1., maxrₑ=400., numrₑ=100, kwargs...)
    """
    Compute the line profile for a given metric using a 
    thin disk and lamppost corona.
    """

    # Setting the inner radius to the ISCO if the entered value is <0
    if minrₑ < 0.
        minrₑ = Gradus.isco(m)
    end

    # Disk
    d = ThinDisc(minrₑ, Inf)

    # Setting up the model and emissivity profile
    model = LampPostModel(h = height)
    profile = emissivity_profile(m, d, model)

    # Computing the line profile
    _, flux = lineprofile(m, x, d, profile; verbose=true, bins=bins, 
            method=TransferFunctionMethod(), minrₑ=minrₑ, maxrₑ=maxrₑ, numrₑ=200)

    return flux
end

function RenderImage(m, x; λ_max, imageSize=(40,30), kwargs...)
    """
    Render a redshift image of the disk
    """

    # Disk
    d = ThinDisc(0.0, 15.0)

    # Redshift point function
    redshift = ConstPointFunctions.redshift(m, x)
    redshiftGeometry = redshift ∘ ConstPointFunctions.filter_intersected()

    # Rendering the image
    α, β, image = rendergeodesics(
            m, x, d, λ_max, pf = redshiftGeometry,
            # image parameters
            image_width = imageSize[1], image_height = imageSize[2],
            αlims = (-20, 20), βlims = (-15, 15), verbose = true)

    return Dict("α" => α, "β" => β, "image" => image)
end

function JohannsenParamVar(setup, bins, func; kwargs...)
    """
    Either render a redshift image or compute a line profile for the Johannsen metric 
    with given parameters and a lamppost corona at a given height above the accretion disk.
    """

    # Position of the observer
    x = SVector(0.0, 10000.0, deg2rad(setup["θ"]), 0.0)

    # Instantiating the metric
    m = JohannsenMetric(setup["M"], setup["a"], setup["α13"], setup["α22"], setup["α52"], setup["ϵ3"])

    # Calling the selected function
    output = func(m, x; λ_max=2x[2], height = setup["h"], bins=bins, kwargs...)

    return output
end

function ParamLoop(parameter, values, setupDict; bins=collect(range(0.1, 1.5, 180)), 
                   line=true, render=true, ParamVar=JohannsenParamVar, kwargs...)
    """
    Loop through a list of values for a given parameter and either render a redshift image,
    compute a line profile, or both.
    """

    # Setup arrays for storage
    results = Dict("flux" => [], "α" => [], "β" => [], "image" => [], "config" => [])

    # Looping through the list of values
    for i in values
        try
            # Setting the value
            setupDict[parameter] = i
            config = generateConfig([parameter], setupDict)

            # Computing line profile
            if line
                flux = ParamVar(setupDict, bins, ComputeLineProfile; kwargs...)
                append!(results["flux"], [flux])
            end

            # rendering redshift image
            if render
                result = ParamVar(setupDict, bins, RenderImage; kwargs...)

                for key in ["α", "β", "image"]
                    append!(results[key], [result[key]])
                end
            end

            # Storing the configuration
            append!(results["config"], [config])
        
        catch err
            println("Value of $i failed to compute")
            println(err)
        end
    end

    # Plotting the results
    plotting(parameter, results, bins; line, render)
end

function Plotting(parameter, results, bins; line, render)
    """
    Plot the results from ParamLoop
    """

    # Pulling the configs used
    configs = results["config"]

    # Plotting the redshift images
    if render
        # Plotting the heatmaps as subplots
        display(heatmap(results["α"], results["β"], results["image"], aspect_ratio = 1; 
                layout = length(configs), title=[i for j in 1:1, i in configs], titleloc=:center, cmap=:redsblues))

        # Plotting the heatmaps individually
        for i in eachindex(configs)
            display(heatmap(results["α"][i], results["β"][i], results["image"][i], aspect_ratio = :equal, title=configs[i]; clims=(0., 1.5), cmap=:redsblues))
            savefig("testimages/$(configs[i]).png")
        end
    end

    # Plotting the line profiles
    if line
        # Plotting the line profile
        plot(xlabel="g", ylabel="Flux (Arbitrary)")

        for i in eachindex(configs)
            plot!(bins, results["flux"][i]; label=configs[i], palette=:tab10, lw=3)
        end

        display(plot!(title="Variations of $parameter"))
    end
end