using Gradus
using Plots
using ColorSchemes
using SpectralFitting

include("LampPostModelFit.jl")

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

function paramVar(setup, bins, func; kwargs...)

    # Position of the observer
    x = SVector(0.0, 10000.0, deg2rad(setup["θ"]), 0.0)
    λ_max = 2x[2]  # Max affine time ~2x[2]

    # Instantiating the metric
    m = JohannsenMetric(setup["M"], setup["a"], setup["α13"], setup["α22"], setup["α52"], setup["ϵ3"])

    output = func(m, x; λ_max=λ_max, height = setup["h"], bins=bins, kwargs...)

    return output

end

function paramLoop(parameter, values, setupDict; bins=collect(range(0.1, 1.5, 180)), 
                   line=true, render=true, kwargs...)

    # Setup arrays for storage
    results = Dict("flux" => [],
                   "α" => [],
                   "β" => [],
                   "image" => [],
                   "config" => [])

    for i in values
        try
            setupDict[parameter] = i

            config = generateConfig([parameter], setupDict)
            if line
                flux = paramVar(setupDict, bins, ComputeLineProfile; kwargs...)
                append!(results["flux"], [flux])
            end

            if render
                result = paramVar(setupDict, bins, RenderImage; kwargs...)

                for key in ["α", "β", "image"]
                    append!(results[key], [result[key]])
                end
            end

            append!(results["config"], [config])

        catch err
            println("Value of $i failed to compute")
            println(err)
        end
    end

    plotting(parameter, results, bins; line, render)
end

function plotting(parameter, results, bins; line, render)

    configs = results["config"]

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

    if line
        # Plotting the line profile
        plot(xlabel="g", ylabel="Flux (Arbitrary)")

        for i in eachindex(configs)
            plot!(bins, results["flux"][i]; label=configs[i], palette=:tab10, lw=3)
        end

        display(plot!(title="Variations of $parameter"))
    end
end

# =======================================================================
# Setup
# =======================================================================

# Metric Defaults
const M = 1.0     # Mass
const a = 0.998   # Spin

# Perturbation Defaults
const α13 = 0.
const α22 = 0.
const α52 = 0.
const ϵ3  = 0.

# BH Defaults
const θ = 60. # Inclination, degrees
const h = 10. # Corona height

# Dictionary for easy access and modification of the model parameters
defaultSetupDict = Dict((["θ", θ], ["α13", α13], ["M", M], ["α22", α22], ["ϵ3", ϵ3], ["a", a], ["h", h], ["α52", α52]))

# =======================================================================
# Function calls
# =======================================================================

# Looping through the values for each variable and computing the line profile
# paramLoop("ϵ3", 8:1:10, copy(defaultSetupDict); line=true, render=true)

bins = collect(range(0.2, 1.5, 100))
flux = paramVar(defaultSetupDict, bins, ComputeLineProfile; 
                render = false, minrₑ = -1., maxrₑ = 400., 
                numrₑ = 100)

plot(data, markersize=3)
display(plot!(xlabel="Energy", ylabel="Flux (arb. units)"))

data = InjectiveData(bins, flux, name="Data")

flux += rand(-0.001:1e-7:0.001, length(flux))
flux[flux.<0] .= 0

model = LampPostJohannsen()

prob = FittingProblem(model => data)

println("+ Fitting...")
result = SpectralFitting.fit(prob, LevenbergMarquadt(); autodiff = :finite, maxIter = 15, verbose=true)

plot!(result)
